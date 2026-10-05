"""
make_mock_hfss_frequencies.py
Fabricate a multi-frequency HFSS data tree from the single real 500 GHz
dataset, to exercise the frequency-keyed lookup.

For each id, frequency f (position i in --freqs) and Ephi in {0,1}, copies
  <src>/<id>_<source-freq>GHz_Ephi=<n>/{far_field,waveguide}.csv
to
  <dst>/waveguides/<id>_<f>GHz_Ephi=<n>/
applying two per-frequency signatures so that which dataset the simulation
actually sampled is observable in the output:

  1. transmittance  - waveguide.csv OutgoingPower scaled by --scales[i], so the
     crack transmittance differs per frequency;
  2. exit direction - for ODD i, far-field amplitudes are zeroed wherever
     Theta > 90 deg. Theta is measured from the gap axis, so a transmitted
     photon from those sets always exits with k_z >= 0, while even-i sets keep
     both signs. This catches a mix-up that scaling alone cannot see, because
     the direction CDF is normalised and a global scale is invisible to it;
  3. sidecar        - each mock frequency gets <id>_<f>GHz.dataset.json, a copy
     of the source's with the frequency, propagating_count, file checksums and
     provenance (mock_of, transmittance_scale, far_field_theta_above_90_zeroed)
     rewritten.

The Freq column of both files is rewritten to "<f>GHz"; the C++ loader checks
it against the directory name (BBR009). Every other column is copied verbatim
(line-by-line text, not a parse/reformat round trip), so the entry whose scale
is 1.0 reproduces the real data exactly.

Writing into --src, or anywhere under the data tree that holds it (the
parent of --src, i.e. data/), is refused: mock data must never enter data/.

Usage:
  conda run -n bbrsim python validation/Scripts/make_mock_hfss_frequencies.py \\
      --src data/waveguides --dst build/mock_hfss \\
      --ids InfParallelPlate_crack1Rohan InfParallelPlate_crack2 \\
      --source-freq 500 --freqs 50 150 500 1500 5000 --scales 0.2 0.4 1.0 0.6 0.8
"""
import argparse
import copy
import os
import sys
import time

from bbrsim import sidecar

# Column indices (both files start with Freq; headers are fixed, see the CSVs).
FF_THETA = 5                    # far_field.csv: Theta [deg]
FF_AMPS = (6, 7, 8, 9)          # rEphi_real, rEphi_imag, rEtheta_real, rEtheta_imag
WG_OUTPOWER = 4                 # waveguide.csv: OutgoingPower

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True,
                help="real waveguides dir (holds <id>_<source-freq>GHz_Ephi=N)")
ap.add_argument("--dst", required=True,
                help="mock data ROOT; datasets are written under <dst>/waveguides/")
ap.add_argument("--ids", nargs="+", required=True)
ap.add_argument("--source-freq", default="500", help="frequency token of the real data")
ap.add_argument("--freqs", nargs="+", default=["50", "150", "500", "1500", "5000"])
ap.add_argument("--scales", nargs="+", type=float, default=[0.2, 0.4, 1.0, 0.6, 0.8])
args = ap.parse_args()

src = os.path.realpath(args.src)
dst = os.path.realpath(args.dst)
# The datasets land in <dst>/waveguides, so guard that directory, and keep the
# whole data tree holding --src (data/, which cmake --install copies) clean.
out_root = os.path.realpath(os.path.join(dst, "waveguides"))
data_root = os.path.dirname(src)


def inside(path, root):
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


if inside(out_root, src) or inside(dst, data_root):
    print(f"refusing: --dst {dst} writes into the real data tree {data_root}; "
          "mock data must not enter it")
    sys.exit(2)
if len(args.freqs) != len(args.scales):
    print(f"--freqs ({len(args.freqs)}) and --scales ({len(args.scales)}) "
          "must have the same length")
    sys.exit(2)
src_sidecars = {}
for id_ in args.ids:
    try:
        src_sidecars[id_] = sidecar.load(src, f"{id_}_{args.source_freq}GHz")
    except ValueError as e:
        print(f"source sidecar: {e}")
        sys.exit(2)


def transform(src_path, dst_path, freq_token, scale, zero_back_hemisphere, is_far_field):
    """Copy a CSV line by line, rewriting only Freq, OutgoingPower and amplitudes."""
    with open(src_path) as fin, open(dst_path, "w") as fout:
        fout.write(fin.readline())                     # header, unchanged
        for line in fin:
            if not line.strip():
                continue
            v = line.rstrip("\n").split(",")
            v[0] = f"{freq_token}GHz"
            if is_far_field:
                if zero_back_hemisphere and float(v[FF_THETA]) > 90.0:
                    for i in FF_AMPS:
                        v[i] = "0.0"
            elif scale != 1.0:
                v[WG_OUTPOWER] = repr(float(v[WG_OUTPOWER]) * scale)
            fout.write(",".join(v) + "\n")


t0 = time.time()
n_files = 0
for id_ in args.ids:
    for i, (ftoken, scale) in enumerate(zip(args.freqs, args.scales)):
        zero_back = (i % 2 == 1)
        for ephi in (0, 1):
            sdir = os.path.join(src, f"{id_}_{args.source_freq}GHz_Ephi={ephi}")
            ddir = os.path.join(out_root, f"{id_}_{ftoken}GHz_Ephi={ephi}")
            if not os.path.isdir(sdir):
                print(f"missing source dataset {sdir}")
                sys.exit(2)
            if os.path.exists(ddir) and os.path.samefile(ddir, sdir):
                print(f"refusing: {ddir} is the source dataset {sdir}")
                sys.exit(2)
            os.makedirs(ddir, exist_ok=True)
            for name, is_ff in (("far_field.csv", True), ("waveguide.csv", False)):
                transform(os.path.join(sdir, name), os.path.join(ddir, name),
                          ftoken, scale, zero_back, is_ff)
                n_files += 1
        stem = f"{id_}_{ftoken}GHz"
        sc = copy.deepcopy(src_sidecars[id_])
        f_ghz = float(ftoken)
        sc["frequency_ghz"] = f_ghz
        sc["frequency_label"] = f"{ftoken}GHz"
        sc["provenance"] = {**sc["provenance"], "producer": "validation/Scripts/make_mock_hfss_frequencies.py",
                            "mock_of": f"{id_}_{args.source_freq}GHz", "transmittance_scale": scale,
                            "far_field_theta_above_90_zeroed": zero_back}
        sc["modes"]["propagating_count"] = sidecar.modes_for(sc["exit_field"]["cross_section"], f_ghz)["propagating_count"]
        sc["files"] = sidecar.file_entries(out_root, stem)
        sidecar.write(out_root, stem, sc)
        print(f"  {id_}_{ftoken}GHz  scale={scale:<4g} "
              f"far-field={'k_z>=0 only' if zero_back else 'full'}")
print(f"wrote {n_files} files to {out_root} in {time.time() - t0:.1f} s")
