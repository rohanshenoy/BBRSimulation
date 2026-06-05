"""
check_reflectance.py
Analyse build/test_output.csv produced by reflectance.mac and compare
observed reflectance R_obs against the Drude model prediction R_theory.

Usage:
    conda run -n bbrsim python scripts/check_reflectance.py [--csv path]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

# ── physical constants (SI) ───────────────────────────────────────────────────
eps0     = 8.8541878128e-12   # F/m
m_e      = 9.109e-31          # kg
n_e      = 8.49e28            # m^-3
e_C      = 1.602e-19          # C
c_light  = 2.998e8            # m/s
sigma_RT = 5.96e7             # S/m (universal Cu at 273 K)
h_eVs    = 4.13566769692e-15  # eV·s

def drude_sigma_dc(RRR, T_K):
    """DC conductivity via Matthiessen's rule (phonons frozen below 50 K)."""
    sigma_imp = RRR * sigma_RT
    if T_K >= 50.0:
        sigma_ph = sigma_RT * 273.0 / T_K
        return 1.0 / (1.0/sigma_imp + 1.0/sigma_ph)
    return sigma_imp

def drude_R(energy_eV, RRR, T_K):
    """Normal-incidence Fresnel reflectance from the Drude model."""
    sigma_dc = drude_sigma_dc(RRR, T_K)
    tau      = sigma_dc * m_e / (n_e * e_C**2)
    nu       = energy_eV / h_eVs          # Hz
    omega    = 2.0 * np.pi * nu
    ot       = omega * tau
    sigma_r  = sigma_dc / (1.0 + ot**2)
    ratio    = sigma_r / (eps0 * omega)
    root     = np.sqrt(1.0 + ratio**2)
    k_       = (omega / c_light) * np.sqrt(0.5 * (root + 1.0))
    kap      = (omega / c_light) * np.sqrt(0.5 * (root - 1.0))
    n_re     = c_light * k_  / omega
    n_im     = c_light * kap / omega
    R = ((n_re - 1.0)**2 + n_im**2) / ((n_re + 1.0)**2 + n_im**2)
    return float(np.clip(R, 0.0, 1.0))

# ── parse args ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--csv", default=os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "build", "test_output.csv"))
parser.add_argument("--RRR",  type=int,   default=100)
parser.add_argument("--T_K",  type=float, default=4.0)
parser.add_argument("--freq", type=float, default=500.0,
                    help="Gun frequency [GHz] (default 500)")
args = parser.parse_args()

if not os.path.exists(args.csv):
    print(f"ERROR: {args.csv} not found.  Run reflectance.mac first.", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(args.csv)

# ── identify Cu boundary hits ─────────────────────────────────────────────────
# With Drude model, Cu materials are named "Cu_RRR{N}_T{T}K".
cu_mask = df["mat_post"].str.startswith("Cu_RRR", na=False)
df_cu   = df[cu_mask].copy()

if df_cu.empty:
    print("No Cu boundary steps found in CSV.  Did the gun hit the Cu slab?")
    sys.exit(1)

# ── classify events: absorbed vs reflected ────────────────────────────────────
# For the reflectance.mac geometry (normal incidence, no cracks in path):
#   Absorbed:  Cu boundary is the ONLY boundary step for that track (max n_reflect == 1).
#   Reflected: Cu boundary + world-exit step → max n_reflect == 2.
max_nr  = df.groupby("event_id")["n_reflect"].max()
cu_evts = set(df_cu["event_id"].unique())

n_hit      = len(cu_evts)
n_absorbed = sum(1 for eid in cu_evts if max_nr[eid] == 1)
n_reflect  = n_hit - n_absorbed

R_obs  = n_reflect / n_hit if n_hit > 0 else float("nan")
D_obs  = 1.0 - R_obs

# ── theory prediction ─────────────────────────────────────────────────────────
energy_eV = args.freq * 1e9 * h_eVs     # GHz → Hz → eV
R_theory  = drude_R(energy_eV, args.RRR, args.T_K)
D_theory  = 1.0 - R_theory

sigma_dc  = drude_sigma_dc(args.RRR, args.T_K)
tau_ps    = sigma_dc * m_e / (n_e * e_C**2) * 1e12
f_break   = 1.0 / (2.0 * np.pi * (tau_ps * 1e-12)) / 1e9
wt        = 2.0 * np.pi * args.freq * 1e9 * (tau_ps * 1e-12)

# ── report ────────────────────────────────────────────────────────────────────
print(f"\n{'─'*60}")
print(f"  BBRsim reflectance check — Cu_RRR{args.RRR}_T{int(args.T_K)}K")
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

# Statistical uncertainty on R_obs
sigma_R = np.sqrt(R_obs * (1.0 - R_obs) / n_hit) if n_hit > 0 else float("nan")
pull    = (R_obs - R_theory) / sigma_R if sigma_R > 0 else float("nan")
print(f"  σ(R)     = {sigma_R:.2e}   pull = {pull:+.2f} σ")

tol = 5.0
print()
if abs(pull) < tol:
    print(f"  PASS  (|pull| = {abs(pull):.2f} σ  <  {tol} σ)")
else:
    print(f"  FAIL  (|pull| = {abs(pull):.2f} σ  ≥  {tol} σ)")
    sys.exit(1)
print(f"{'─'*60}\n")
