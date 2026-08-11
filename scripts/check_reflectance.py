"""
check_reflectance.py
Analyse build/output/bbr.root produced by reflectance.mac and compare
observed reflectance R_obs against the Drude model prediction R_theory.

Usage:
    conda run -n bbrsim python scripts/check_reflectance.py [--root path]

`--csv` is accepted as a deprecated alias for `--root` (the input has been a
ROOT file since the Phase-A output migration; the flag name is historical).
"""

import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings
from bbrsim import physics, select

# ── parse args ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--root", "--csv", dest="root", default=os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "build", "output", "bbr.root"),
    help="Path to bbr.root (--csv is a deprecated alias)")
parser.add_argument("--RRR",  type=int,   default=100)
parser.add_argument("--T_K",  type=float, default=4.0)
parser.add_argument("--freq", type=float, default=500.0,
                    help="Gun frequency [GHz] (default 500)")
args = parser.parse_args()

if not os.path.exists(args.root):
    print(f"ERROR: {args.root} not found.  Run reflectance.mac first.", file=sys.stderr)
    sys.exit(1)

df = load_crossings(args.root)

if select.cu_boundary(df).empty:
    print("No Cu boundary steps found.  Did the gun hit the Cu slab?")
    sys.exit(1)

n_hit, n_absorbed = select.cu_absorption_stats(df)
n_reflect = n_hit - n_absorbed

R_obs  = n_reflect / n_hit if n_hit > 0 else float("nan")
D_obs  = 1.0 - R_obs

# ── theory prediction (from bbrsim.physics) ──────────────────────────────────
freq_Hz   = args.freq * 1e9
R_theory  = float(physics.drude_reflectance(freq_Hz, args.RRR, args.T_K))
D_theory  = 1.0 - R_theory

sigma_dc  = physics.sigma_dc(args.RRR, args.T_K)
tau_s     = physics.drude_tau(args.RRR, args.T_K)
tau_ps    = tau_s * 1e12
f_break   = 1.0 / (2.0 * np.pi * tau_s) / 1e9
wt        = 2.0 * np.pi * freq_Hz * tau_s

# ── report ────────────────────────────────────────────────────────────────────
print(f"\n{'─'*60}")
print(f"  BBRsim reflectance check — Cu_RRR{args.RRR}_T{args.T_K:g}K")
print(f"{'─'*60}")
print(f"  σ_DC    = {sigma_dc:.3e} S/m")
print(f"  τ       = {tau_ps:.3f} ps")
print(f"  f_break = {f_break:.1f} GHz")
print(f"  ωτ at {args.freq:.0f} GHz = {wt:.2f}  "
      f"{'(Drude needed: ωτ > 1)' if wt > 1 else '(H-R valid: ωτ < 1)'}")
print()
print(f"  Events total     : {len(df['event_id'].unique())}")
print(f"  Cu boundary hits : {n_hit}")
print(f"  Absorbed         : {n_absorbed}")
print(f"  Reflected        : {n_reflect}")
print()
print(f"  R_obs    = {R_obs:.6f}   (D_obs    = {D_obs:.4e})")
print(f"  R_theory = {R_theory:.6f}   (D_theory = {D_theory:.4e})")

# Statistical test on the absorbed COUNT (Poisson). With the full Drude model
# D_theory ~ 5e-5, so the expected number of absorptions in a 10k-event run is
# O(1) and the binomial pull on R_obs degenerates (σ→0 when R_obs = 1).
lam = n_hit * D_theory
if lam > 0:
    pull = (n_absorbed - lam) / np.sqrt(lam)
else:
    pull = float("nan")
print(f"  expected absorbed (λ = N·D_theory) = {lam:.2f}   observed = {n_absorbed}")
print(f"  Poisson pull = {pull:+.2f} σ")

tol = 5.0
print()
if abs(pull) < tol:
    print(f"  PASS  (|pull| = {abs(pull):.2f} σ  <  {tol} σ)")
else:
    print(f"  FAIL  (|pull| = {abs(pull):.2f} σ  ≥  {tol} σ)")
    sys.exit(1)
print(f"{'─'*60}\n")
