"""
make_mock_schema2.py
Relabel the real 500 GHz crack datasets as schema 2.0, for the schema-2 reader's end-to-end case
(OR-028, OR-034, OR-038). Each legacy dataset
  <src>/<id>_<freq>GHz_Ephi=<n>/{far_field,waveguide}.csv
becomes
  <dst>/waveguides/<id>_<freq>GHz_Ephi=<n>/{far_field,waveguide}.csv
  <dst>/waveguides/<id>_<freq>GHz.dataset.json          (schema 2.0)
with the 2.0 waveguide columns: the header is relabelled (the legacy files store the six field
components interleaved, Ex_re, Ex_im, Ey_re, ..., HF-038, so the field values stay as they are, position
by position), IngoingPower is unchanged, and InterceptFactor is appended, c_k of the 15-key grid
(|cos theta|, cell_90 at 90 deg). The far-field files are copied byte for byte.

A basis is dark where its key's source ratio t = OutgoingPower / IngoingPower satisfies |t| <= 1e-12
(the three theta = 0 keys per polarization); a ratio below -1e-12 is refused. The legacy theta = 0
powers are tiny but nonzero (1e-44), which OR-038 forbids on a dark basis, so OutgoingPower is written
as exactly 0 on every row of a dark basis; lit bases keep the source OutgoingPower text. The stored
fields are all preserved, W = 1, and the sidecar is built (bbrsim.sidecar.build_schema2_from_csvs) from the
finalized files, so its checksums are theirs. Provenance: mock, of <legacy stem>, scale_rule_holds false
(the legacy fields were never scaled to the exit flux), material vacuum.

Writing into --src, or anywhere under the data tree that holds it, is refused: mock data must never
enter data/.

Usage:
  conda run -n bbrsim python validation/Scripts/make_mock_schema2.py --src data/waveguides \\
      --dst build/mock_v2 --ids InfParallelPlate_crack1Rohan InfParallelPlate_crack2
"""
import argparse
import os
import shutil
import sys
import time

from bbrsim import sidecar

DARK_REL = 1e-12
OUT_COL, IN_COL, PHI_COL, THETA_COL = 4, 5, 2, 3

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True, help="real waveguides dir (holds <id>_<freq>GHz_Ephi=N and sidecars)")
ap.add_argument("--dst", required=True, help="mock data ROOT; datasets are written under <dst>/waveguides/")
ap.add_argument("--ids", nargs="+", required=True)
ap.add_argument("--source-freq", default="500", help="frequency token of the real data")
args = ap.parse_args()

src, dst = os.path.realpath(args.src), os.path.realpath(args.dst)
out_root = os.path.realpath(os.path.join(dst, "waveguides"))
data_root = os.path.dirname(src)


def inside(path, root):
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


if inside(out_root, src) or inside(dst, data_root):
    print(f"refusing: --dst {dst} writes into the real data tree {data_root}; mock data must not enter it")
    sys.exit(2)


def key_of(v):
    return (round(float(v[PHI_COL]), 2), round(float(v[THETA_COL]), 2))


def relabel(src_path, dst_path, thetas):
    """Copy waveguide.csv line by line into the 2.0 layout; returns {key: dark}."""
    dark = {}
    with open(src_path) as fin, open(dst_path, "w") as fout:
        fin.readline()
        fout.write(",".join(sidecar.WG_COLUMNS_2) + "\n")
        for line in fin:
            if not line.strip():
                continue
            v = line.rstrip("\n").split(",")
            if len(v) != 15:
                print(f"{src_path}: a row has {len(v)} fields, not 15")
                sys.exit(2)
            k = key_of(v)
            if k not in dark:
                t = float(v[OUT_COL]) / float(v[IN_COL])
                if t < -DARK_REL:
                    print(f"{src_path}: key {k} has a negative transmittance {t!r}")
                    sys.exit(2)
                dark[k] = abs(t) <= DARK_REL
            if dark[k]:
                v[OUT_COL] = "0.0"
            ck = sidecar.intercept_factor(k[1], thetas)
            fout.write(",".join(v) + "," + repr(ck) + "\n")
    return dark


t0 = time.time()
os.makedirs(out_root, exist_ok=True)
for id_ in args.ids:
    stem = f"{id_}_{args.source_freq}GHz"
    try:
        s = sidecar.load(src, stem)
    except ValueError as e:
        print(f"source sidecar: {e}")
        sys.exit(2)
    thetas = s["excitation"]["incident_theta_deg"]
    n_dark = 0
    for ephi in (0, 1):
        sdir = os.path.join(src, f"{stem}_Ephi={ephi}")
        ddir = os.path.join(out_root, f"{stem}_Ephi={ephi}")
        if not os.path.isdir(sdir):
            print(f"missing source dataset {sdir}")
            sys.exit(2)
        if os.path.exists(ddir) and os.path.samefile(ddir, sdir):
            print(f"refusing: {ddir} is the source dataset {sdir}")
            sys.exit(2)
        os.makedirs(ddir, exist_ok=True)
        shutil.copyfile(os.path.join(sdir, "far_field.csv"), os.path.join(ddir, "far_field.csv"))
        dark = relabel(os.path.join(sdir, "waveguide.csv"), os.path.join(ddir, "waveguide.csv"), thetas)
        n_dark += sum(dark.values())
    ex, geo = s["exit_field"], s["geometry"]
    sc = sidecar.build_schema2_from_csvs(
        out_root, id_, f"{args.source_freq}GHz", section=ex["cross_section"], extent_mm=geo["extent_mm"],
        provenance={"producer": "validation/Scripts/make_mock_schema2.py", "mock": True, "of": stem,
                    "scale_rule_holds": False, "material": "vacuum"},
        exit_origin_mm=s["frames"]["exit_cs"]["origin_mm_global"],
        plane_wave_origin_mm=s["excitation"]["origin_mm_global"], bounding_box_mm=geo["bounding_box_mm"],
        resolution_mm=ex["resolution_mm"], outside_points=ex["outside_points"], rim_points=ex["rim_points"],
        bounds_method=ex["bounds_method"], rotational=s["symmetry"]["rotational"], walls=s["boundaries"]["walls"])
    sidecar.write(out_root, f"{stem}", sc)
    sidecar.check_full(out_root, id_, stem, float(args.source_freq))
    print(f"  {stem}: schema 2.0, {n_dark} dark (key, polarization) pairs set to OutgoingPower 0")
print(f"wrote {len(args.ids)} schema-2 datasets to {out_root} in {time.time() - t0:.1f} s")
