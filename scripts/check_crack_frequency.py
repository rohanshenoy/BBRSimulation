"""
check_crack_frequency.py
Validate the frequency-keyed HFSS lookup on the output of crack_frequency.mac:
16 runs, one ROOT file each (output/bbr_freq_rNN.root), taken against the MOCK
five-frequency tree written by scripts/make_mock_hfss_frequencies.py.

What makes a wrong implementation fail here:
  * recording a frequency but sampling another dataset - the mock scales the
    transmittance per frequency (weakest pair separated by 7 sigma at 4000
    photons), so T_obs pins which dataset was really used;
  * mixing far_field.csv and waveguide.csv across frequencies, which the
    transmittance alone cannot see - the 150 and 1500 GHz mock sets have their
    far field truncated to Theta <= 90 deg, so their transmitted photons must
    all exit with k_z >= 0 while the others use both signs;
  * nearest-in-LINEAR frequency instead of log - the eight probes sit 1 % on
    either side of the four log midpoints, where the two rules disagree;
  * a first-selection-sticks cache, or per-photon staleness - the broadband
    Planck run is checked photon by photon;
  * a missing -1 sentinel - every non-crack crossing must carry it.

Clamping (runs 13-14) changes only the warning, never the selected dataset, so
it is checked by the runner counting BBR008 lines, not here. The tie rule is
unobservable from a macro; scripts/check_physics.py unit-tests it instead.

Usage:
    conda run -n bbrsim python scripts/check_crack_frequency.py [output_dir] [--data-dir MOCK_ROOT]
    defaults: build/output, and <output_dir>/../../mock_hfss
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import hfss, select
from bbrsim.io import load_many

GRID = [50.0, 150.0, 500.0, 1500.0, 5000.0]
MIDS = [np.sqrt(a * b) for a, b in zip(GRID[:-1], GRID[1:])]
CRACK1 = "InfParallelPlate_crack1Rohan"
# Mock sets whose far field was truncated to Theta <= 90 deg (odd positions in
# the generator's --freqs list): transmitted photons exit with k_z >= 0 only.
KZ_POSITIVE_SETS = {150.0, 1500.0}
# Planck bin floors: half the expectation from the spec's emitter Monte Carlo
# (about 200 / 1500 / 8600 / 12900 / 800 entries at 100k events).
PLANCK_MIN = {50.0: 100, 150.0: 700, 500.0: 4000, 1500.0: 6000, 5000.0: 400}
SIGMA_T = 4.0

# Run order of crack_frequency.mac: (label, nu_GHz or None for Planck,
# expected selection or None, events).
RUNS = [("near 50 GHz", 50.05, 50.0, 4000),
        ("on-grid 150 GHz", 150.0, 150.0, 4000),
        ("on-grid 500 GHz", 500.0, 500.0, 4000),
        ("on-grid 1500 GHz", 1500.0, 1500.0, 4000),
        ("near 5000 GHz", 4995.0, 5000.0, 4000)]
for m, lo, hi in zip(MIDS, GRID[:-1], GRID[1:]):
    RUNS.append((f"{m:.1f} x0.99", m * 0.99, lo, 4000))
    RUNS.append((f"{m:.1f} x1.01", m * 1.01, hi, 4000))
RUNS.append(("20 GHz (clamp low)", 20.0, 50.0, 4000))
RUNS.append(("10 THz (clamp high)", 10000.0, 5000.0, 4000))
RUNS.append(("Planck 20 K", None, None, 100000))

ap = argparse.ArgumentParser()
ap.add_argument("output_dir", nargs="?", default="build/output")
ap.add_argument("--data-dir", default=None,
                help="mock data ROOT (contains waveguides/); "
                     "default <output_dir>/../../mock_hfss")
args = ap.parse_args()
mock_root = args.data_dir or os.path.normpath(
    os.path.join(args.output_dir, "..", "..", "mock_hfss"))
wg_dir = os.path.join(mock_root, "waveguides")

files = sorted(glob.glob(os.path.join(args.output_dir, "bbr_freq_r*.root")))
if not files:
    print(f"no bbr_freq_r*.root under {args.output_dir}")
    print("RESULT: FAIL")
    sys.exit(1)

df = load_many(files)
entries = select.crack_crossings(df)
DIFFRACTED = ["BBRDiffractionTransmit", "BBRDiffractionReflect"]
decided = entries[entries["status"].isin(DIFFRACTED)]

results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))


print(f"files: {len(files)}   crossings: {len(df)}   crack entries: {len(decided)}")
print(f"mock tree: {wg_dir}")

# --- 1. discovery mirrors the C++ for every crack the run actually touched ---
ids = sorted({v.split(":")[0] for v in decided["vol_post"].dropna().unique()})
grids = {}
for cid in ids:
    ent = hfss.discover_frequencies(cid, wg_dir)
    grids[cid] = ent
    check(f"discovery {cid}", [f for f, _ in ent] == GRID, f"{[f for f, _ in ent]}")

# --- 2. column present, grid value on crack entries, -1 everywhere else ------
has_col = "hfss_freq_GHz" in df.columns
check("hfss_freq_GHz column present", has_col)
if not has_col:
    for n, ok, d in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {n:<46} {d}")
    print("RESULT: FAIL")
    sys.exit(1)
others = df.drop(decided.index)
n_bad_sentinel = int((others["hfss_freq_GHz"] != -1.0).sum())
check("sentinel -1 off crack entries", n_bad_sentinel == 0,
      f"{n_bad_sentinel} of {len(others)} rows carry a frequency")
check("grid value on every crack entry", decided["hfss_freq_GHz"].isin(GRID).all())

# --- 3. per-photon selection equals the Python mirror ------------------------
nu = hfss.photon_frequency_GHz(decided["energy_eV"].to_numpy(float))
expect = np.array([hfss.select_frequency(grids[v.split(":")[0]], n)[0]
                   for v, n in zip(decided["vol_post"], nu)])
n_mismatch = int((decided["hfss_freq_GHz"].to_numpy(float) != expect).sum())
check("selection == python mirror, every photon", n_mismatch == 0,
      f"{n_mismatch} mismatches of {len(decided)}")

# --- per-run behaviour --------------------------------------------------------
crack1 = decided[decided["vol_post"].str.contains("crack1", na=False)]
print(f"\n{'run':>3} {'config':<20} {'N':>6} {'sel':>6} {'T_obs':>7} {'T_pred':>7} "
      f"{'z':>6}  exit k_z")
for r, (label, nu_r, f_exp, n_ev) in enumerate(RUNS):
    rows = crack1[crack1["run_id"] == r]
    if nu_r is None:                                   # broadband Planck run
        counts = rows["hfss_freq_GHz"].value_counts()
        for f, floor in PLANCK_MIN.items():
            check(f"run{r:02d} planck bin {f:g} GHz >= {floor}",
                  counts.get(f, 0) >= floor, f"{int(counts.get(f, 0))} entries")
        print(f"{r:>3} {label:<20} {len(rows):>6}   bins "
              f"{ {int(k): int(v) for k, v in sorted(counts.items())} }")
        continue

    N = len(rows)
    check(f"run{r:02d} one entry per event",
          N == n_ev and (rows.groupby("event_id").size() == 1).all(),
          f"{N} entries for {n_ev} events")
    if N == 0:
        continue
    sel = rows["hfss_freq_GHz"].unique()
    check(f"run{r:02d} selects {f_exp:g} GHz", len(sel) == 1 and sel[0] == f_exp,
          f"{sorted(sel)}")
    f_sel = float(sel[0])

    # Transmittance against the mock dataset THIS run selected.
    stem = dict(grids[CRACK1])[f_sel]
    ds = hfss.load_dataset(stem, wg_dir)
    inc = hfss.fold_incidence(rows[["px_pre", "py_pre", "pz_pre"]].to_numpy(float)[0])
    key = hfss.nearest_key(ds, inc.phi_deg, inc.theta_deg)
    T_pred, _ = hfss.random_polarization_mixture(ds[key])
    tx = rows[rows["status"] == "BBRDiffractionTransmit"]
    T_obs = len(tx) / N
    sigma = np.sqrt(T_pred * (1 - T_pred) / N)
    z = (T_obs - T_pred) / sigma
    check(f"run{r:02d} T_obs matches the selected mock set", abs(z) < SIGMA_T,
          f"{T_obs:.4f} vs {T_pred:.4f} ({z:+.2f} sigma)")
    check(f"run{r:02d} transmitted + reflected == entries",
          len(tx) + len(rows[rows["status"] == "BBRDiffractionReflect"]) == N)

    # Far-field signature: which CSV the direction sampler really read.
    pz = tx["pz_post"].to_numpy(float)
    if f_sel in KZ_POSITIVE_SETS:
        check(f"run{r:02d} far field truncated (k_z >= 0)", (pz >= -1e-12).all(),
              f"min pz {pz.min():+.4f}")
        sig = "+ only"
    else:
        check(f"run{r:02d} far field full (both k_z signs)",
              (pz > 0).any() and (pz < 0).any(),
              f"min {pz.min():+.3f} max {pz.max():+.3f}")
        sig = "both"
    if f_sel == 500.0:                                  # the unscaled real data
        check(f"run{r:02d} reproduces the real 50 %",
              abs(T_obs - 0.5) < SIGMA_T * np.sqrt(0.25 / N), f"{T_obs:.4f}")

    print(f"{r:>3} {label:<20} {N:>6} {f_sel:>6g} {T_obs:>7.4f} {T_pred:>7.4f} "
          f"{z:>+6.2f}  {sig}")

print()
n_fail = sum(1 for _, ok, _ in results if not ok)
for name, ok, detail in results:
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<46} {detail}")
print(f"\nchecks: {len(results)}   failed: {n_fail}")
print(f"RESULT: {'PASS' if n_fail == 0 else 'FAIL'}")
sys.exit(0 if n_fail == 0 else 1)
