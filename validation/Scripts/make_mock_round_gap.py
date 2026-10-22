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
  * a 2 um exit lattice (1961 points per key at the default R) keeping the twenty lattice
    points on r = R, the solutions of a^2 + b^2 = 25^2 such as (30, 40) um
    (rim_points included), which BBRsim clamps 1 nm inside the disc at load; the lattice
    spans +-R, so --radius must be a multiple of 2 um (exit 2 otherwise);
  * a far-field Gaussian lobe (sigma 25 deg) around (Theta, Phi) = (90, Phi0), with
    Phi0 = -30 deg at IWaveTheta = 135 and 0 elsewhere.
Refuses to write inside the data tree that holds --real, and exits 2 when --real
holds no *.dataset.json. A re-run repoints crack links that name another --real.

--schema 2 writes the same dataset natively in schema 2.0 (OR-028, OR-034, OR-038): the 16 waveguide
columns with the six field components in their stored (interleaved) positions Ex_re, Ex_im, Ey_re,
Ey_im, Ez_re, Ez_im, InterceptFactor c_k of the declared grid (|cos theta|, cell_90 at 90 deg), and
OutgoingPower = T * IngoingPower * c_k so the loaded T equals the table above; the dark keys
(T = 0) hold exactly 0. Provenance carries material "vacuum". The cracks are not linked: the schema-2
cracks come from make_mock_schema2.py written into the same --dst, and --real is then optional, used
only for the refusal to write into a data tree.
"""
import argparse
import math
import os
import shutil
import sys

from bbrsim import sidecar

ap = argparse.ArgumentParser()
ap.add_argument("--real", help="real waveguides dir (crack datasets with sidecars); required for schema 1")
ap.add_argument("--schema", type=int, choices=(1, 2), default=1, help="sidecar schema of the mock (default 1)")
ap.add_argument("--dst", required=True, help="mock data ROOT; the tree is written under <dst>/waveguides/")
ap.add_argument("--id", default="RoundGap_r50um")
ap.add_argument("--freq", default="2000")
ap.add_argument("--radius", type=float, default=5e-5, help="HFSS radius [m]")
ap.add_argument("--length", type=float, default=4e-4, help="gap length [m]")
args = ap.parse_args()

if args.real is None and args.schema == 1:
    ap.error("--real is required for --schema 1")
dst = os.path.realpath(args.dst)
out = os.path.join(dst, "waveguides")
if args.real is not None:
    real = os.path.realpath(args.real)
    data_root = os.path.dirname(real)
    if dst == data_root or dst.startswith(data_root.rstrip(os.sep) + os.sep):
        print(f"refusing: --dst {dst} lies inside the real data tree {data_root}")
        sys.exit(2)
    REAL_SIDECARS = sorted(n for n in os.listdir(real) if n.endswith(".dataset.json")) if os.path.isdir(real) else []
    if not REAL_SIDECARS and args.schema == 1:
        print(f"error: --real {real} holds no *.dataset.json; pass the real waveguides directory (data/waveguides)")
        sys.exit(2)
if args.schema == 2 or args.real is None:
    REAL_SIDECARS = []                           # schema 2 links no cracks (make_mock_schema2.py writes them)

R, L = args.radius, args.length
STEP = 2e-6                                                # exit lattice step [m]
N_HALF = round(R / STEP) if math.isfinite(R) and R > 0 else 0
if N_HALF < 1 or abs(N_HALF * STEP - R) > 1e-12 * R:
    print(f"error: --radius {R!r} m is not a positive multiple of the {STEP!r} m exit-lattice step; "
          "the lattice must reach r = R so the rim points exist")
    sys.exit(2)
MU0 = 1.25663706212e-6
INGOING = math.pi * R * R / (2 * sidecar.C * MU0)          # Ei = 1 V/m over the entrance disc
PHIS, THETAS = [0.0, 45.0, 90.0], [0.0, 45.0, 90.0, 135.0, 180.0]
T = {180.0: (0.8, 0.6), 135.0: (0.5, 0.3), 90.0: (0.1, 0.05), 45.0: (0.0, 0.0), 0.0: (0.0, 0.0)}
FF_THETA = [15.0 * i for i in range(13)]                   # 0 .. 180
FF_PHI = [-90.0 + 15.0 * i for i in range(12)]             # -90 .. 75
LATTICE = [round(2 * k * 1e-6, 12) for k in range(-N_HALF, N_HALF + 1)]   # 2 um step, +-R
EXIT = [(y, z) for y in LATTICE for z in LATTICE if y * y + z * z <= R * R * (1 + sidecar.SECTION_REL)]
label = f"{args.freq}GHz"
stem = f"{args.id}_{label}"


CK = {t: sidecar.intercept_factor(t, THETAS) for t in THETAS}        # schema 2: c_k per key
WG_HEADER = sidecar.WG_COLUMNS_2 if args.schema == 2 else sidecar.WG_COLUMNS


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
        fh.write(",".join(WG_HEADER) + "\n")
        for phi in PHIS:
            for theta in THETAS:
                if max_t(phi, theta) > 1.0:
                    sys.exit(f"table error: max T > 1 at ({phi}, {theta})")
                if args.schema == 2:
                    outp = T[theta][e] * INGOING * CK[theta]       # T = OutgoingPower / (IngoingPower * c_k)
                    for y, z in EXIT:
                        ex, ey = fields(e, phi, y, z)
                        fh.write(f"{label},{e},{phi},{theta},{outp!r},{INGOING!r},0.0,{y!r},{z!r},"
                                 f"{ex!r},0.0,{ey!r},0.0,0.0,0.0,{CK[theta]!r}\n")
                    continue
                outp = T[theta][e] * INGOING
                for y, z in EXIT:
                    ex, ey = fields(e, phi, y, z)
                    fh.write(f"{label},{e},{phi},{theta},{outp!r},{INGOING!r},0.0,{y!r},{z!r},"
                             f"{ex!r},{ey!r},0.0,0.0,0.0,0.0\n")

for name in REAL_SIDECARS:
    real_stem = name[: -len(".dataset.json")]
    shutil.copy2(os.path.join(real, name), os.path.join(out, name))
    for e in (0, 1):
        link, target = os.path.join(out, f"{real_stem}_Ephi={e}"), os.path.join(real, f"{real_stem}_Ephi={e}")
        if os.path.islink(link):
            if os.readlink(link) == target:
                continue
            os.remove(link)                                # a link from an earlier --real: repoint it
        elif os.path.lexists(link):
            print(f"error: {link} exists and is not a symlink; remove it or choose another --dst")
            sys.exit(2)
        os.symlink(target, link)

build = sidecar.build_schema2_from_csvs if args.schema == 2 else sidecar.build_from_csvs
mock_provenance = {"producer": "validation/Scripts/make_mock_round_gap.py", "hand_written": False, "mock": True,
                   "of": "cylindrical2 (radius 50 um, length 0.4 mm, 2000 GHz), synthetic tables"}
if args.schema == 2:
    mock_provenance = {**mock_provenance, "scale_rule_holds": False, "material": "vacuum"}
sc = build(
    out, args.id, label, section={"shape": "disc", "radius_m": R},
    extent_mm={"p": L * 1e3, "l": 2 * R * 1e3, "g": 2 * R * 1e3},
    provenance=mock_provenance,
    exit_origin_mm=[0.0, 0.0, L * 1e3], plane_wave_origin_mm=[0.0, 0.0, 0.0],
    bounding_box_mm=[-R * 1e3, -R * 1e3, 0.0, R * 1e3, R * 1e3, L * 1e3],
    resolution_mm=[0.0, 0.002, 0.002], outside_points="omitted", rim_points="included", rotational=True)
sidecar.write(out, stem, sc)
sidecar.check_full(out, args.id, stem, float(args.freq))
if args.schema == 2:
    print(f"wrote {stem} in schema 2.0 ({len(EXIT)} exit points per key) into {out}")
else:
    print(f"wrote {stem} ({len(EXIT)} exit points per key) and linked the real cracks into {out}")
