#!/usr/bin/env python3
# scripts/check_reflectance.py
# Analytical Hagen-Rubens verification for OFHC Cu at RRR=100, 4 K.
# Run from repo root: conda run -n bbrsim python scripts/check_reflectance.py

import numpy as np

eps0   = 8.8541878128e-12   # F/m
sigma  = 5.96e9             # S/m  (RRR=100 × σ_RT(Cu) = 5.96e7)
h_eVs  = 4.13566769692e-15  # eV·s

print(f"{'Frequency':>10}  {'R':>12}  {'A=1-R':>10}")
for nu_GHz, label in [(50,"50 GHz"),(500,"500 GHz"),(5000,"5 THz"),(20000,"20 THz")]:
    nu    = nu_GHz * 1e9
    omega = 2 * np.pi * nu
    R     = 1.0 - 2.0 * np.sqrt(2.0 * eps0 * omega / sigma)
    R     = float(np.clip(R, 0.0, 1.0))
    print(f"{label:>10}  {R:.8f}  {1-R:.4e}")

nu500  = 500e9
R500   = 1.0 - 2.0 * np.sqrt(2.0 * eps0 * 2*np.pi*nu500 / sigma)
assert 0.999 < R500 < 1.000, f"R(500 GHz) = {R500:.8f} out of [0.999, 1.000]"

R50 = 1.0 - 2.0 * np.sqrt(2.0 * eps0 * 2*np.pi*50e9 / sigma)
assert R50 > R500, "R should increase toward lower frequencies"

print("PASS: Hagen-Rubens formula verified for OFHC Cu RRR=100")
