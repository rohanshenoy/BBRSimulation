"""
check_invariants.py
Validate the two invariants any BBRsim output must satisfy, one printed section
each: no photon travels inside a metal, and the abspoints.term_status labels
agree with the crossings ntuple. PASS only if both hold.

No photons in metal. An optical photon never propagates INSIDE a REFLECTIVITY
metal: BBSimOpBoundaryProcess either reflects it or kills it at the surface, so
no boundary-crossing row may start in a metal (mat_pre must never be a
Cu_RRR* / BBR_Perfect* material). The file must also hold at least one
crossing, so an empty output cannot pass; --allow-no-crossings waives only this
requirement, for a fixture that crosses no boundary by design
(Validation_WorldExit.mac). Runs with no metal contact at all (e.g.
Validation_CrackTransmit.mac at normal incidence) are still meaningful: a
photon that skated into a crack wall and was reflected into the copper would
appear here as a row starting in Cu.

Termination labels. The file must hold at least one abspoints row, and
  1. no termination row is labelled "unknown";
  2. every photon that left the world (term_vol == "none") is labelled
     "WorldExit";
  3. every row labelled as an absorption (BBRAbsorb, Absorption, Detection,
     BulkAbsorption) terminated inside a volume (term_vol != "none");
  4. the number of BBRAbsorb termination rows equals the number of BBRAbsorb
     crossings rows (one boundary absorption kills exactly one photon).

Breaks this catches:
  - a wrong reflection normal on a daughter->mother crossing (photon reflected
    INTO the copper instead of back into the vacuum);
  - primaries generated inside a metal volume (emitter overlapping a wall);
  - reading the boundary process's per-thread status on a world-exit step,
    where Geant4 never invoked it, so the label is stale from an earlier step
    or track (StepTooSmall, or even BBRAbsorb for a photon that was never
    absorbed), plus "unknown" before the first boundary crossing per thread.

Usage:
    conda run -n bbrsim python validation/check_invariants.py [path/to/bbr.root] [--allow-no-crossings]
"""
import argparse
import sys

from bbrsim.io import load

METAL_PREFIXES = ("Cu_RRR", "BBR_Perfect")
ABSORB = {"BBRAbsorb", "Absorption", "Detection", "BulkAbsorption"}

ap = argparse.ArgumentParser(description="No photons in metal, and correct term_status labels.")
ap.add_argument("path", nargs="?", default="output/bbr.root")
ap.add_argument("--allow-no-crossings", action="store_true",
                help="let the no-photons-in-metal check pass on a file with no "
                     "crossings (a fixture that crosses no boundary by design)")
args = ap.parse_args()
PATH = args.path
cr, ab = load(PATH)


def is_metal(series):
    m = series.str.startswith(METAL_PREFIXES[0], na=False)
    for p in METAL_PREFIXES[1:]:
        m |= series.str.startswith(p, na=False)
    return m


# ---- no photons in metal -----------------------------------------------------
print("== no photons in metal")
in_metal = cr[is_metal(cr["mat_pre"])]
metal_hits = cr[is_metal(cr["mat_post"])]

print(f"crossings rows            : {len(cr)}")
print(f"rows entering a metal     : {len(metal_hits)}")
print(f"rows STARTING in a metal  : {len(in_metal)}   (must be 0)")
if len(in_metal):
    print("  examples (event_id, vol_pre -> vol_post, status):")
    for _, r in in_metal.head(5).iterrows():
        print(f"    {r['event_id']:>7}  {r['vol_pre']} -> {r['vol_post']}  {r['status']}")

passed_metal = len(in_metal) == 0 and (len(cr) > 0 or args.allow_no_crossings)
if len(cr) == 0:
    print("no crossings at all — allowed by --allow-no-crossings" if args.allow_no_crossings
          else "no crossings at all — empty output cannot pass")
print(f"no photons in metal : {'PASS' if passed_metal else 'FAIL'}")

# ---- termination labels ------------------------------------------------------
print("== termination labels")
if len(ab) == 0:
    print("no abspoints rows — nothing to check")
    print("termination labels  : FAIL")
    print("RESULT: FAIL (" + ("" if passed_metal else "no photons in metal, ") + "termination labels)")
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
print(f"termination labels  : {'PASS' if passed else 'FAIL'}")

failed = [s for s, ok in (("no photons in metal", passed_metal), ("termination labels", passed)) if not ok]
passed = not failed
print("RESULT: " + ("FAIL (" + ", ".join(failed) + ")" if failed else "PASS"))
sys.exit(0 if passed else 1)
