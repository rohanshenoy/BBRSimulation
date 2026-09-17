"""
check_term_status.py
Validate the abspoints.term_status labels against the crossings ntuple.

Invariants (any BBRsim output):
  1. no termination row is labelled "unknown";
  2. every photon that left the world (term_vol == "none") is labelled
     "WorldExit";
  3. every row labelled as an absorption (BBRAbsorb, Absorption, Detection,
     BulkAbsorption) terminated inside a volume (term_vol != "none");
  4. the number of BBRAbsorb termination rows equals the number of BBRAbsorb
     crossings rows (one boundary absorption kills exactly one photon).

Break this catches: reading the boundary process's per-thread status on a
world-exit step, where Geant4 never invoked it, so the label is stale from an
earlier step or track (StepTooSmall, or even BBRAbsorb for a photon that was
never absorbed), plus "unknown" before the first boundary crossing per thread.

Usage:
    conda run -n bbrsim python scripts/check_term_status.py [path/to/bbr.root]
"""
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load

ABSORB = {"BBRAbsorb", "Absorption", "Detection", "BulkAbsorption"}

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"
cr, ab = load(PATH)

if len(ab) == 0:
    print("no abspoints rows — nothing to check")
    print("RESULT: FAIL")
    sys.exit(1)

n_unknown   = int((ab["term_status"] == "unknown").sum())
left_world  = ab[ab["term_vol"] == "none"]
n_left_bad  = int((left_world["term_status"] != "WorldExit").sum())
absorbed    = ab[ab["term_status"].isin(ABSORB)]
n_abs_none  = int((absorbed["term_vol"] == "none").sum())
n_abs_term  = int((ab["term_status"] == "BBRAbsorb").sum())
n_abs_cross = int((cr["status"] == "BBRAbsorb").sum())

print(f"abspoints rows                          : {len(ab)}")
print("term_status counts                      :", ab["term_status"].value_counts().to_dict())
print(f"1. 'unknown' labels                     : {n_unknown}   (must be 0)")
print(f"2. world exits not labelled WorldExit   : {n_left_bad} of {len(left_world)}   (must be 0)")
print(f"3. absorption labels with term_vol none : {n_abs_none}   (must be 0)")
print(f"4. BBRAbsorb terminations / crossings   : {n_abs_term} / {n_abs_cross}   (must match)")

passed = n_unknown == 0 and n_left_bad == 0 and n_abs_none == 0 and n_abs_term == n_abs_cross
print(f"RESULT: {'PASS' if passed else 'FAIL'}")
sys.exit(0 if passed else 1)
