"""
check_cu_absorptance.py
Parse BBRSim stdout for [BBR] reflectance lines; compare A_obs to
Hagen-Rubens Planck-weighted theory. PASS if 0.3 < A_obs/A_theory < 3.0.

Usage:
    ./BBRSim test.mac >bbrsim_out.txt 2>&1
    conda run -n bbrsim python scripts/check_cu_absorptance.py bbrsim_out.txt [--sigma S]

Options:
  --sigma S   Effective conductivity in S/m used to compute A_theory.
              Defaults to 5.96e9 (OFHC_Cu RRR=100 theoretical).
              Use 1.786e8 for OF_Cu, 3.383e8 for HP_Cu.
"""

import sys
import re
import argparse
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="?", default="bbrsim_stdout.txt")
parser.add_argument("--sigma", type=float, default=5.96e9,
                    help="sigma_eff in S/m for Hagen-Rubens theory (default: 5.96e9 = OFHC_Cu)")
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

eps0  = 8.854e-12   # F/m
sigma = args.sigma
hbar  = 1.0546e-34  # J·s
k_B   = 1.3806e-23  # J/K
T     = 4.0         # K

def A_HR(omega):
    return 2.0 * np.sqrt(2.0 * eps0 * omega / sigma)

def planck_weight(omega):
    x = hbar * omega / (k_B * T)
    return np.where(x < 500, omega**2 / (np.expm1(x)), 0.0)

omega_lo = 2 * np.pi * 1e9
omega_hi = 2 * np.pi * 20e12

num, _ = integrate.quad(lambda w: A_HR(w) * planck_weight(w), omega_lo, omega_hi)
den, _ = integrate.quad(lambda w:            planck_weight(w), omega_lo, omega_hi)
A_theory = num / den

ratio = A_obs / A_theory
lo, hi = 0.3, 3.0
passed = lo < ratio < hi

print(f"sigma_eff used             = {sigma:.4e} S/m")
print(f"A_theory (Planck-weighted) = {A_theory:.6e}")
print(f"A_obs                      = {A_obs:.6e}")
print(f"ratio A_obs/A_theory       = {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")

sys.exit(0 if passed else 1)
