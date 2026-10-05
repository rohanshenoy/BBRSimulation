"""
make_mock_round_gap.py
A mock HFSS dataset for the straight round gap (Jason Wang's cylindrical2:
radius 50 um, length 0.4 mm, 2000 GHz), in the canonical pose and in sidecar
schema 1.0, for validation/G4Macros/Validation_RoundGap.mac. Writes
  <dst>/waveguides/RoundGap_r50um_2000GHz_Ephi={0,1}/{far_field,waveguide}.csv
  <dst>/waveguides/RoundGap_r50um_2000GHz.dataset.json
and links the real crack datasets of --real beside them (with copies of their
sidecars): BBRsim validates every placed crack, and the test world always
places crack1 and crack2.

The tables are synthetic, chosen so the validator can see each mapping:
  * T per incidence key differs between Ephi=0 and Ephi=1;
  * polarization-specific radial exit profiles, Ephi=0 E = (1 - (r/R)^2, 0, 0) and
    Ephi=1 E = (0.3 (1 - (r/R)^2) at IWavePhi = 45 else 0, (r/R)^2, 0), so rho != 0
    at the diagonal keys and the cross term is active there;
  * a 2 um exit lattice (1961 points per key) that keeps the twenty lattice
    points on r = R, the solutions of a^2 + b^2 = 25^2 such as (30, 40) um
    (rim_points included), so the Geant4 radius needs its 1 um margin;
  * a far-field Gaussian lobe (sigma 25 deg) around (Theta, Phi) = (90, Phi0), with
    Phi0 = -30 deg at IWaveTheta = 135 and 0 elsewhere.
Refuses to write inside the data tree that holds --real.
"""
import argparse
import math
import os
import shutil
import sys

from bbrsim import sidecar

ap = argparse.ArgumentParser()
ap.add_argument("--real", required=True, help="real waveguides dir (crack datasets with sidecars)")
ap.add_argument("--dst", required=True, help="mock data ROOT; the tree is written under <dst>/waveguides/")
ap.add_argument("--id", default="RoundGap_r50um")
ap.add_argument("--freq", default="2000")
ap.add_argument("--radius", type=float, default=5e-5, help="HFSS radius [m]")
ap.add_argument("--length", type=float, default=4e-4, help="gap length [m]")
args = ap.parse_args()

real, dst = os.path.realpath(args.real), os.path.realpath(args.dst)
out, data_root = os.path.join(dst, "waveguides"), os.path.dirname(real)
if dst == data_root or dst.startswith(data_root.rstrip(os.sep) + os.sep):
    print(f"refusing: --dst {dst} lies inside the real data tree {data_root}")
    sys.exit(2)

R, L = args.radius, args.length
MU0 = 1.25663706212e-6
INGOING = math.pi * R * R / (2 * sidecar.C * MU0)          # Ei = 1 V/m over the entrance disc
PHIS, THETAS = [0.0, 45.0, 90.0], [0.0, 45.0, 90.0, 135.0, 180.0]
T = {180.0: (0.8, 0.6), 135.0: (0.5, 0.3), 90.0: (0.1, 0.05), 45.0: (0.0, 0.0), 0.0: (0.0, 0.0)}
FF_THETA = [15.0 * i for i in range(13)]                   # 0 .. 180
FF_PHI = [-90.0 + 15.0 * i for i in range(12)]             # -90 .. 75
LATTICE = [round(k * 1e-6, 12) for k in range(-50, 51, 2)]
EXIT = [(y, z) for y in LATTICE for z in LATTICE if y * y + z * z <= R * R * (1 + sidecar.SECTION_REL)]
label = f"{args.freq}GHz"
stem = f"{args.id}_{label}"


def fields(e, phi, y, z):
    s = (y * y + z * z) / (R * R)
    if e == 0:
        return (1.0 - s, 0.0)
    return (0.3 * (1.0 - s) if phi == 45.0 else 0.0, s)


def max_t(phi, theta):
    """Top eigenvalue of [[T0, c], [c, T1]]: must not exceed 1 (no load-time normalization)."""
    t0, t1 = T[theta]
    a = [fields(0, phi, y, z) for y, z in EXIT]
    b = [fields(1, phi, y, z) for y, z in EXIT]
    x = sum(p[0] * q[0] + p[1] * q[1] for p, q in zip(a, b))
    p0, p1 = sum(p[0] ** 2 + p[1] ** 2 for p in a), sum(q[0] ** 2 + q[1] ** 2 for q in b)
    c = math.sqrt(t0 * t1) * x / math.sqrt(p0 * p1)
    return 0.5 * (t0 + t1) + math.hypot(0.5 * (t0 - t1), c)


os.makedirs(out, exist_ok=True)
for e in (0, 1):
    d = os.path.join(out, f"{stem}_Ephi={e}")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "far_field.csv"), "w") as fh:
        fh.write(",".join(sidecar.FF_COLUMNS) + "\n")
        for phi in PHIS:
            for theta in THETAS:
                phi0 = -30.0 if theta == 135.0 else 0.0
                for t in FF_THETA:
                    for p in FF_PHI:
                        a = math.exp(-((t - 90.0) ** 2 + (p - phi0) ** 2) / (2 * 25.0 ** 2))
                        r_ephi, r_etheta = (0.0, a) if e == 0 else (a, 0.0)
                        fh.write(f"{label},{e},{phi},{theta},{p},{t},{r_ephi!r},0.0,{r_etheta!r},0.0\n")
    with open(os.path.join(d, "waveguide.csv"), "w") as fh:
        fh.write(",".join(sidecar.WG_COLUMNS) + "\n")
        for phi in PHIS:
            for theta in THETAS:
                if max_t(phi, theta) > 1.0:
                    sys.exit(f"table error: max T > 1 at ({phi}, {theta})")
                outp = T[theta][e] * INGOING
                for y, z in EXIT:
                    ex, ey = fields(e, phi, y, z)
                    fh.write(f"{label},{e},{phi},{theta},{outp!r},{INGOING!r},0.0,{y!r},{z!r},"
                             f"{ex!r},{ey!r},0.0,0.0,0.0,0.0\n")

for name in sorted(os.listdir(real)):
    if not name.endswith(".dataset.json"):
        continue
    real_stem = name[: -len(".dataset.json")]
    shutil.copy2(os.path.join(real, name), os.path.join(out, name))
    for e in (0, 1):
        link = os.path.join(out, f"{real_stem}_Ephi={e}")
        if not os.path.lexists(link):
            os.symlink(os.path.join(real, f"{real_stem}_Ephi={e}"), link)

sc = sidecar.build_from_csvs(
    out, args.id, label, section={"shape": "disc", "radius_m": R},
    extent_mm={"p": L * 1e3, "l": 2 * R * 1e3, "g": 2 * R * 1e3},
    provenance={"producer": "validation/Scripts/make_mock_round_gap.py", "hand_written": False, "mock": True,
                "of": "cylindrical2 (radius 50 um, length 0.4 mm, 2000 GHz), synthetic tables"},
    exit_origin_mm=[0.0, 0.0, L * 1e3], plane_wave_origin_mm=[0.0, 0.0, 0.0],
    bounding_box_mm=[-R * 1e3, -R * 1e3, 0.0, R * 1e3, R * 1e3, L * 1e3],
    resolution_mm=[0.0, 0.002, 0.002], outside_points="omitted", rim_points="included", rotational=True)
sidecar.write(out, stem, sc)
sidecar.check_full(out, args.id, stem, float(args.freq))
print(f"wrote {stem} ({len(EXIT)} exit points per key) and linked the real cracks into {out}")
