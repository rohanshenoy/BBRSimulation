"""
check_cu_absorptance.py
Parse BBRSim stdout for [BBR] reflectance lines; compare A_obs to
full Drude Planck-weighted theory. PASS if 0.3 < A_obs/A_theory < 3.0.

Usage:
    ./BBRSim test.mac >bbrsim_out.txt 2>&1
    conda run -n bbrsim python scripts/check_cu_absorptance.py bbrsim_out.txt [--rrr N]

Options:
  --rrr N    Residual Resistance Ratio of the Cu used in the simulation.
             Defaults to 100 (OFHC_Cu, Cu_RRR100_T4K).
             Use 3 for OF_Cu (Cu_RRR3_T4K), 6 for HP_Cu (Cu_RRR6_T4K).
             σ_DC = RRR × σ_RT is derived internally (σ_RT = 5.96e7 S/m).
"""

import sys
import re
import argparse
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="?", default="bbrsim_stdout.txt")
parser.add_argument("--rrr", type=int, default=100,
                    help="RRR of Cu material (default 100 = OFHC_Cu)")
args = parser.parse_args()

pattern = re.compile(
    r"\[BBR\] reflectance mat=(\S+) N=(\d+) A_obs=([\d.e+-]+) R_theory=([\d.e+-]+)"
)
best = None
with open(args.file) as f:
    for line in f:
        m = pattern.search(line)
        if m:
            mat, n, a_obs = m.group(1), int(m.group(2)), float(m.group(3))
            if best is None or n > best[1]:
                best = (mat, n, a_obs)

if best is None:
    print("ERROR: no [BBR] reflectance lines found in", args.file)
    sys.exit(1)

mat_name, N, A_obs = best
print(f"[BBR] log  mat={mat_name}  N={N}  A_obs={A_obs:.6e}")

# Physical constants (SI)
sigma_RT = 5.96e7        # S/m, universal for Cu at 273 K
m_e      = 9.109e-31     # kg
n_e      = 8.49e28       # m^-3
e_C      = 1.602e-19     # C
eps0     = 8.854e-12     # F/m
c        = 2.998e8       # m/s
hbar     = 1.0546e-34    # J·s
k_B      = 1.3806e-23    # J/K
T_K      = 4.0           # K (cryogenic baseline)

# DC conductivity: at 4K impurity dominates, σ_DC = RRR × σ_RT
sigma_DC = args.rrr * sigma_RT
tau      = sigma_DC * m_e / (n_e * e_C * e_C)

def drude_R(omega):
    """Normal-incidence reflectance from full Drude model (Griffiths §9.4)."""
    ot      = omega * tau
    sigma_r = sigma_DC / (1.0 + ot * ot)
    ratio   = sigma_r / (eps0 * omega)
    root    = np.sqrt(1.0 + ratio * ratio)
    k_wav   = (omega / c) * np.sqrt(0.5 * (root + 1.0))
    kappa   = (omega / c) * np.sqrt(0.5 * (root - 1.0))
    n_re    = c * k_wav / omega
    n_im    = c * kappa  / omega
    return np.clip(((n_re - 1.0)**2 + n_im**2) / ((n_re + 1.0)**2 + n_im**2), 0.0, 1.0)

def planck_weight(omega):
    x = hbar * omega / (k_B * T_K)
    return np.where(x < 500, omega**2 / (np.expm1(x)), 0.0)

omega_lo = 2 * np.pi * 1e9     # 1 GHz lower cutoff for integral
omega_hi = 2 * np.pi * 20e12   # 20 THz upper cutoff

num, _ = integrate.quad(lambda w: (1.0 - drude_R(w)) * planck_weight(w), omega_lo, omega_hi)
den, _ = integrate.quad(lambda w:                       planck_weight(w), omega_lo, omega_hi)
A_theory = num / den

ratio = A_obs / A_theory
lo, hi = 0.3, 3.0
passed = lo < ratio < hi

print(f"RRR                        = {args.rrr}")
print(f"sigma_DC                   = {sigma_DC:.4e} S/m")
print(f"tau                        = {tau * 1e12:.3f} ps")
print(f"A_theory (Planck+Drude)    = {A_theory:.6e}")
print(f"A_obs                      = {A_obs:.6e}")
print(f"ratio A_obs/A_theory       = {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")

sys.exit(0 if passed else 1)
