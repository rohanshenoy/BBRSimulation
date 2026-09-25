"""
check_crack_transmittance.py
Validate the HFSS crack transmittance at normal incidence on the output of
Validation_CrackTransmit.mac (40 000 photons, 500 GHz, into crack1, seeds 2024 7).

Expected physics (hand-derived): a sub-cutoff parallel-plate gap is a perfect
polarization filter — the TEM component transmits (T=1), the orthogonal
component is cut off (T=0) — so for unpolarized light T = 0.50 exactly.

Two breaks this catches:
  - raw HFSS power ratios > 1 (port-normalization artefact, 1.0545 for crack1)
    used without renormalization: T_obs is ~2.2 points high, a ~9 sigma miss;
  - far-field samples lying exactly in the exit-face plane (HFSS Phi = -90 deg
    grid edge): a transmitted photon leaves with zero normal component.

PASS if |T_obs - 0.50| < 3 sigma_binomial and no transmitted photon has
|dir . n| < 1e-6, and every event entered crack1 exactly once.

Usage:
    conda run -n bbrsim python validation/check_crack_transmittance.py [path/to/bbr.root]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings
from bbrsim import select

T_IDEAL = 0.50
PATH = sys.argv[1] if len(sys.argv) > 1 else "output/bbr.root"

df = load_crossings(PATH)
n_events = df["event_id"].nunique()
entries = select.crack_crossings(df)
entries = entries[entries["vol_post"].str.contains("crack1", na=False)]
one_entry_each = (entries.groupby("event_id").size() == 1).all() and len(entries) == n_events

tx = entries[entries["status"] == "BBRDiffractionTransmit"]
N, k = len(entries), len(tx)
T_obs = k / N if N else float("nan")
sigma = np.sqrt(T_IDEAL * (1 - T_IDEAL) / N) if N else float("nan")
z = (T_obs - T_IDEAL) / sigma if N else float("nan")

# Normal component of the exit direction: crack1's exit face normal is +x.
n_tangential = int((np.abs(tx["px_post"].values) < 1e-6).sum())

print(f"events                     : {n_events}")
print(f"crack1 entries             : {N}   (one per event: {'yes' if one_entry_each else 'NO'})")
print(f"transmitted                : {k}")
print(f"T_obs                      : {T_obs:.5f} +/- {sigma:.5f}")
print(f"T_ideal                    : {T_IDEAL:.2f}   (z = {z:+.2f} sigma, threshold 3)")
print(f"tangential exits |px|<1e-6 : {n_tangential}   (must be 0)")

# Frequency-keyed lookup (2026-09-22): each crack entry records the grid
# frequency the wrapper selected. The real data hold only 500 GHz, so every
# entry must read 500 and every other crossing the -1 sentinel.
has_col = "hfss_freq_GHz" in df.columns
decided = df[(df["mat_post"] == "vacuum_wg")
             & df["status"].isin(["BBRDiffractionTransmit", "BBRDiffractionReflect"])] \
    if has_col else df.iloc[0:0]
freq_ok = bool(has_col and len(decided)
               and (decided["hfss_freq_GHz"] == 500.0).all()
               and (df.drop(decided.index)["hfss_freq_GHz"] == -1.0).all())
print(f"hfss_freq_GHz              : {'present' if has_col else 'MISSING'}; "
      f"500 on {len(decided)} crack entries, -1 elsewhere: "
      f"{'yes' if freq_ok else 'NO'}")

passed = bool(one_entry_each and abs(z) < 3.0 and n_tangential == 0 and freq_ok)
print(f"RESULT: {'PASS' if passed else 'FAIL'}")
sys.exit(0 if passed else 1)
