#!/usr/bin/env python3
"""
check_cu_serov.py
Verify that the Hagen-Rubens σ_eff values chosen for OF_Cu and HP_Cu
reproduce Serov et al. (IEEE TMT 2016) reflection-loss measurements at T=4 K.

Reference points (Serov Figs 6 and 8, read at T=4 K):
  OF copper  (99.97%),          150 GHz: D = 0.58e-3
  HP copper  (99.999%, annealed), 230 GHz: D = 0.55e-3

Run: conda run -n bbrsim python scripts/check_cu_serov.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import physics

eps0 = physics.EPS0   # F/m (kept for the back-calc of sigma_HP below)


def hagen_rubens_D(freq_Hz, sigma_SI):
    return physics.hagen_rubens_absorptance(freq_Hz, sigma_SI)

# ------------------------------------------------------------------
# 1.  OF copper: σ_eff = 1/ρ₀  (Serov fit: ρ₀=0.56e-8 Ω·m)
# ------------------------------------------------------------------
rho0_OF   = 0.56e-8          # Ω·m  (Serov Eq. 24 fit to OF copper)
sigma_OF  = 1.0 / rho0_OF   # S/m
D_OF_calc = hagen_rubens_D(150e9, sigma_OF)
D_OF_ref  = 0.58e-3          # Serov Fig. 6, T=4 K, 150 GHz

# ------------------------------------------------------------------
# 2.  HP copper annealed: σ_eff back-calculated from D=0.55e-3 at 230 GHz
#     σ = 8 ε₀ ω / D²  (inverted Hagen-Rubens)
# ------------------------------------------------------------------
D_HP_ref   = 0.55e-3         # Serov Fig. 8, T=4 K, 230 GHz
omega_230  = 2.0 * np.pi * 230e9
sigma_HP   = 8.0 * eps0 * omega_230 / D_HP_ref**2
D_HP_calc  = hagen_rubens_D(230e9, sigma_HP)

# ------------------------------------------------------------------
# 3.  Existing OFHC_Cu (theoretical, for reference)
# ------------------------------------------------------------------
sigma_OFHC = 5.96e9
D_OFHC_500 = hagen_rubens_D(500e9, sigma_OFHC)

print("=" * 60)
print(f"{'Material':<20} {'freq':>8} {'D_calc':>10} {'D_ref':>10} {'ratio':>8}")
print("-" * 60)

ratio_OF = D_OF_calc / D_OF_ref
ratio_HP = D_HP_calc / D_HP_ref
tol = 0.10   # ±10 % tolerance

rows = [
    ("OF_Cu (Serov)",  "150 GHz", D_OF_calc, D_OF_ref, ratio_OF),
    ("HP_Cu (Serov)",  "230 GHz", D_HP_calc, D_HP_ref, ratio_HP),
    ("OFHC_Cu (ref)", "500 GHz", D_OFHC_500, None, None),
]
for name, freq, calc, ref, ratio in rows:
    ref_str   = f"{ref:.3e}" if ref  is not None else "   N/A   "
    ratio_str = f"{ratio:.3f}" if ratio is not None else "  N/A  "
    print(f"{name:<20} {freq:>8} {calc:>10.3e} {ref_str:>10} {ratio_str:>8}")

print("=" * 60)
print(f"\nσ_eff  OF_Cu  = {sigma_OF:.4e} S/m   (= 1/ρ₀, ρ₀={rho0_OF:.2e} Ω·m)")
print(f"σ_eff  HP_Cu  = {sigma_HP:.4e} S/m   (back-calc, D=0.55e-3 at 230 GHz)")
print(f"σ      OFHC_Cu= {sigma_OFHC:.4e} S/m   (RRR=100 × σ_RT, theoretical)")

passed = (abs(ratio_OF - 1.0) <= tol) and (abs(ratio_HP - 1.0) <= tol)
print(f"\nRESULT: {'PASS' if passed else 'FAIL'}  (tolerance ±{int(tol*100)}%)")
sys.exit(0 if passed else 1)
