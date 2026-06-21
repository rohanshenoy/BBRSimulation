"""Self-test for analysis/bbrsim/physics.py.

Anchors the shared physics against documented reference numbers (CLAUDE.md,
Serov) so the Python theory cannot silently drift from the C++ BBRMaterials.
PASS only if every anchor is within tolerance.

Run: conda run -n bbrsim python scripts/check_physics.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import physics

ok = True


def check(name, got, expected, rel_tol):
    global ok
    rel = abs(got - expected) / expected
    if rel > rel_tol:
        ok = False
    print(f"  {'PASS' if rel <= rel_tol else 'FAIL'}  {name:<22} "
          f"got={got:.3e}  expected={expected:.3e}  rel={rel:.1%}  tol={rel_tol:.0%}")


print("Drude absorptance at 500 GHz, 4 K (CLAUDE.md reference values):")
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

print()
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
