"""
check_no_photons_in_metal.py
Physics invariant: an optical photon never propagates INSIDE a REFLECTIVITY
metal. BBSimOpBoundaryProcess either reflects it or kills it at the surface,
so no boundary-crossing row may start in a metal (mat_pre must never be a
Cu_RRR* / BBR_Perfect* material).

Breaks this catches:
  - a wrong reflection normal on a daughter->mother crossing (photon reflected
    INTO the copper instead of back into the vacuum);
  - primaries generated inside a metal volume (emitter overlapping a wall).

PASS if zero rows have a metal mat_pre and the file holds at least one
crossing (so an empty output cannot pass). Runs with no metal contact at all
(e.g. crack_transmit.mac at normal incidence) are still meaningful: a photon
that skated into a crack wall and was reflected into the copper would appear
here as a row starting in Cu.

Usage:
    conda run -n bbrsim python scripts/check_no_photons_in_metal.py [path/to/bbr.root]
"""
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings

METAL_PREFIXES = ("Cu_RRR", "BBR_Perfect")

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"
df = load_crossings(PATH)


def is_metal(series):
    m = series.str.startswith(METAL_PREFIXES[0], na=False)
    for p in METAL_PREFIXES[1:]:
        m |= series.str.startswith(p, na=False)
    return m


in_metal = df[is_metal(df["mat_pre"])]
metal_hits = df[is_metal(df["mat_post"])]

print(f"crossings rows            : {len(df)}")
print(f"rows entering a metal     : {len(metal_hits)}")
print(f"rows STARTING in a metal  : {len(in_metal)}   (must be 0)")
if len(in_metal):
    print("  examples (event_id, vol_pre -> vol_post, status):")
    for _, r in in_metal.head(5).iterrows():
        print(f"    {r['event_id']:>7}  {r['vol_pre']} -> {r['vol_post']}  {r['status']}")

passed = len(in_metal) == 0 and len(df) > 0
if len(df) == 0:
    print("no crossings at all — empty output cannot pass")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")
sys.exit(0 if passed else 1)
