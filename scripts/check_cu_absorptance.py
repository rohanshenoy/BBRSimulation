"""
check_cu_absorptance.py
Parse BBRSim output for [BBR] reflectance lines; compare A_obs to
Hagen-Rubens Planck-weighted theory. PASS if 0.3 < A_obs/A_theory < 3.0.

Note: [BBR] lines use G4cout and appear on stdout (not stderr) in Geant4 MT mode.

Usage:
    ./BBRSim test.mac >bbrsim_stdout.txt 2>&1
    conda run -n bbrsim python scripts/check_cu_absorptance.py bbrsim_stdout.txt
"""

import sys
import re
import numpy as np
from scipy import integrate

STDERR_FILE = sys.argv[1] if len(sys.argv) > 1 else "bbrsim_stderr.txt"

# --- Parse [BBR] lines ---
pattern = re.compile(
    r"\[BBR\] reflectance mat=OFHC_Cu N=(\d+) A_obs=([\d.e+-]+) R_theory=([\d.e+-]+)"
)
best = None
with open(STDERR_FILE) as f:
    for line in f:
        m = pattern.search(line)
        if m:
            n, a_obs, r_theory = int(m.group(1)), float(m.group(2)), float(m.group(3))
            if best is None or n > best[0]:
                best = (n, a_obs, r_theory)

if best is None:
    print("ERROR: no [BBR] reflectance lines found in", STDERR_FILE)
    sys.exit(1)

N, A_obs, R_theory = best
print(f"[BBR] log  N={N}  A_obs={A_obs:.6e}  R_theory={R_theory:.6f}")

# --- Compute Planck-weighted <A_theory> at T=4K ---
eps0  = 8.854e-12   # F/m
sigma = 5.96e9      # S/m  (OFHC Cu RRR=100)
hbar  = 1.0546e-34  # J·s
k_B   = 1.3806e-23  # J/K
T     = 4.0         # K

def A_HR(omega):
    return 2.0 * np.sqrt(2.0 * eps0 * omega / sigma)

def planck_weight(omega):
    x = hbar * omega / (k_B * T)
    return np.where(x < 500, omega**2 / (np.expm1(x)), 0.0)

# Integrate from 1 GHz to 20 THz (covers >99.9% of Planck weight at 4K)
omega_lo = 2 * np.pi * 1e9
omega_hi = 2 * np.pi * 20e12

num, _ = integrate.quad(lambda w: A_HR(w) * planck_weight(w), omega_lo, omega_hi)
den, _ = integrate.quad(lambda w:            planck_weight(w), omega_lo, omega_hi)
A_theory_planck = num / den

ratio = A_obs / A_theory_planck
lo, hi = 0.3, 3.0
passed = lo < ratio < hi

print(f"A_theory (Planck-weighted 4K) = {A_theory_planck:.6e}")
print(f"A_obs                         = {A_obs:.6e}")
print(f"ratio A_obs/A_theory          = {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")

if not passed:
    sys.exit(1)
