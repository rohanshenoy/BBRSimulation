"""Self-test for tools/python/bbrsim/physics.py.

Anchors the shared physics against documented reference numbers (user_guide.md,
Serov) so the Python theory cannot silently drift from the C++ BBRMaterials.
PASS only if every anchor is within tolerance.

Run: conda run -n bbrsim python validation/check_physics.py
"""
import os
import sys

import numpy as np

from bbrsim import hfss, paths, physics

ok = True


def check(name, got, expected, rel_tol):
    global ok
    rel = abs(got - expected) / expected
    if rel > rel_tol:
        ok = False
    print(f"  {'PASS' if rel <= rel_tol else 'FAIL'}  {name:<22} "
          f"got={got:.3e}  expected={expected:.3e}  rel={rel:.1%}  tol={rel_tol:.0%}")


def check_bool(name, cond, detail=""):
    global ok
    if not cond:
        ok = False
    print(f"  {'PASS' if cond else 'FAIL'}  {name:<42} {detail}")


print("Drude absorptance at 500 GHz, 4 K (user_guide.md reference values):")
check("RRR=100 D", physics.drude_absorptance(500e9, 100, 4.0), 4.9e-5, 0.15)
check("RRR=3   D", physics.drude_absorptance(500e9,   3, 4.0), 1.0e-3, 0.20)
check("RRR=6   D", physics.drude_absorptance(500e9,   6, 4.0), 6.3e-4, 0.20)

print("Hagen-Rubens vs Serov OF copper (150 GHz, sigma = 1/0.56e-8 S/m):")
check("OF_Cu D", physics.hagen_rubens_absorptance(150e9, 1.0 / 0.56e-8), 0.58e-3, 0.10)

print("Planck photon-number spectrum peak (u = E/kT):")
T = 4.0
u = np.linspace(0.05, 8.0, 4000)
pdf = physics.planck_photon_number_pdf(u * physics.K_EV * T, T)
check("peak u", float(u[int(np.argmax(pdf))]), physics.PLANCK_PEAK_U, 0.02)

print("HFSS frequency selection rule (nearest in log frequency):")
# The C++ (BBRCrackLibrary::Lookup) cannot be exercised for ties or edges from a
# macro, so this is where the rule itself is pinned; the mock-data validator
# then checks that the C++ agrees photon by photon.
GRID = [(50.0, "id_50GHz"), (150.0, "id_150GHz"), (500.0, "id_500GHz"),
        (1500.0, "id_1500GHz"), (5000.0, "id_5000GHz")]
check_bool("single-entry grid never clamps",
           hfss.select_frequency([(500.0, "id_500GHz")], 20.0) == (500.0, "id_500GHz", 0))
check_bool("exact low edge is on-grid",
           hfss.select_frequency(GRID, 50.0) == (50.0, "id_50GHz", 0))
check_bool("exact high edge is on-grid",
           hfss.select_frequency(GRID, 5000.0) == (5000.0, "id_5000GHz", 0))
check_bool("below the grid clamps low",
           hfss.select_frequency(GRID, 20.0) == (50.0, "id_50GHz", -1))
check_bool("above the grid clamps high",
           hfss.select_frequency(GRID, 10000.0) == (5000.0, "id_5000GHz", +1))
for lo, hi in zip([50.0, 150.0, 500.0, 1500.0], [150.0, 500.0, 1500.0, 5000.0]):
    mid = np.sqrt(lo * hi)
    check_bool(f"{mid:8.3f} GHz x 0.99 -> {lo:g}",
               hfss.select_frequency(GRID, mid * 0.99)[0] == lo)
    check_bool(f"{mid:8.3f} GHz x 1.01 -> {hi:g}",
               hfss.select_frequency(GRID, mid * 1.01)[0] == hi)
    check_bool(f"{mid:8.3f} GHz exact tie -> lower {lo:g}",
               hfss.select_frequency(GRID, mid)[0] == lo)
check_bool("photon_frequency_GHz(2.067834e-3 eV) = 500 GHz",
           abs(hfss.photon_frequency_GHz(2.067834e-3) - 500.0) < 1e-3,
           f"got {hfss.photon_frequency_GHz(2.067834e-3):.4f}")

print("bbrsim.paths.data_dir (Python twin of the C++ data default):")
_saved = os.environ.pop("BBRSIMDATA", None)
try:
    os.environ["BBRSIMDATA"] = "/nonexistent/bbrsim-data"
    check_bool("BBRSIMDATA wins", paths.data_dir() == "/nonexistent/bbrsim-data")
    del os.environ["BBRSIMDATA"]
    _d = paths.data_dir()
    # The nearest ancestor of the imported package holding data/waveguides, found
    # independently: a farther tree or another checkout's data must not pass.
    _up, _want = os.path.dirname(os.path.abspath(paths.__file__)), None
    while _want is None:
        if os.path.isdir(os.path.join(_up, "data", "waveguides")):
            _want = os.path.join(_up, "data")
        elif os.path.dirname(_up) == _up:
            break
        _up = os.path.dirname(_up)
    if _want is None and os.path.isdir(os.path.join(sys.prefix, "share", "BBRsim", "data")):
        _want = os.path.join(sys.prefix, "share", "BBRsim", "data")
    check_bool("without BBRSIMDATA: nearest data/waveguides tree",
               _want is not None and _d == _want, f"{_d} (expected {_want})")
finally:
    if _saved is not None:
        os.environ["BBRSIMDATA"] = _saved

print()
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
