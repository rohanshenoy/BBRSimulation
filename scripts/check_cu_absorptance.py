"""
check_cu_absorptance.py
Compare observed Cu absorptance A_obs (from output/bbr.root) to the full-Drude
Planck-weighted theory. A_obs = fraction of Cu boundary crossings that were
absorbed (status BBRAbsorb). RRR/T are inferred from the Cu material name in
the data (override with --rrr / --temp). PASS if 0.3 < A_obs/A_theory < 3.0.

Usage:
    ./BBRSim test.mac          # a Planck run with enough Cu absorptions
    conda run -n bbrsim python scripts/check_cu_absorptance.py [path/to/bbr.root]
                                 [--rrr N] [--temp T]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings
from bbrsim import physics, select

parser = argparse.ArgumentParser()
parser.add_argument("path", nargs="?", default="build/output/bbr.root")
parser.add_argument("--rrr", type=int, default=None,
                    help="Override RRR (default: inferred from Cu material name)")
parser.add_argument("--temp", type=float, default=None,
                    help="Override stage T in K (default: inferred from material name)")
args = parser.parse_args()

if not os.path.exists(args.path):
    print(f"ERROR: {args.path} not found. Run a Planck macro (e.g. test.mac) first.",
          file=sys.stderr)
    sys.exit(1)

df = load_crossings(args.path)
cu = select.cu_boundary(df)
if cu.empty:
    print("ERROR: no Cu boundary crossings found. Did photons reach the Cu slab?",
          file=sys.stderr)
    sys.exit(1)

# Infer RRR/T from the Cu material name (e.g. Cu_RRR100_T4K) unless overridden.
mat_name = cu["mat_post"].iloc[0]
parsed = select.parse_cu_rrr_t(mat_name)
RRR = args.rrr if args.rrr is not None else (parsed[0] if parsed else 100)
T_K = args.temp if args.temp is not None else (parsed[1] if parsed else 4.0)

n_hit, n_absorbed = select.cu_absorption_stats(df)
A_obs = n_absorbed / n_hit if n_hit > 0 else float("nan")
A_theory = physics.planck_weighted_absorptance(RRR, T_K)

ratio = A_obs / A_theory if A_theory > 0 else float("nan")
lo, hi = 0.3, 3.0
passed = lo < ratio < hi

print(f"material           : {mat_name}  (RRR={RRR}, T={T_K:g} K)")
print(f"Cu boundary hits   : {n_hit}")
print(f"absorbed           : {n_absorbed}")
print(f"A_obs              : {A_obs:.6e}")
print(f"A_theory (Planck+Drude) : {A_theory:.6e}")
print(f"ratio A_obs/A_theory    : {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")

sys.exit(0 if passed else 1)
