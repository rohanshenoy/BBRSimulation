#!/usr/bin/env python3
"""
check_cu_serov.py
Verify that the RRR values behind the named Cu aliases (OF_Cu = RRR 3,
HP_Cu = RRR 6; see BBRMaterials::GetCopperByName) reproduce the Serov et al.
(IEEE TMT 2016) reflection-loss measurements at T = 4 K with the full complex
Drude model that BBRsim actually uses.

Reference points (Serov Figs 6 and 8, read at T = 4 K):
  OF copper  (99.97%),           150 GHz: D = 0.58e-3
  HP copper  (99.999%, annealed), 230 GHz: D = 0.55e-3

The expected values are the measured literals; nothing here is derived from
the model under test, so a wrong RRR mapping or a broken Drude formula fails.

Run: conda run -n bbrsim python scripts/check_cu_serov.py
"""
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import physics

T_K = 4.0
TOL = 0.10   # +-10 %

ROWS = [
    # alias, RRR, freq_Hz, D_ref (Serov)
    ("OF_Cu", 3, 150e9, 0.58e-3),
    ("HP_Cu", 6, 230e9, 0.55e-3),
]

print("=" * 64)
print(f"{'Material':<10} {'RRR':>4} {'freq':>9} {'D_Drude':>10} {'D_Serov':>10} {'ratio':>7}")
print("-" * 64)
passed = True
for alias, rrr, f, d_ref in ROWS:
    d_calc = float(physics.drude_absorptance(f, rrr, T_K))
    ratio = d_calc / d_ref
    ok = abs(ratio - 1.0) <= TOL
    passed &= ok
    print(f"{alias:<10} {rrr:>4} {f/1e9:>6.0f} GHz {d_calc:>10.3e} {d_ref:>10.3e} {ratio:>7.3f}  {'ok' if ok else 'OUT OF TOLERANCE'}")

d_ofhc = float(physics.drude_absorptance(500e9, 100, T_K))
print("-" * 64)
print(f"{'OFHC_Cu':<10} {100:>4} {'500 GHz':>9} {d_ofhc:>10.3e} {'(ref only)':>10}")
print("=" * 64)
print(f"\nRESULT: {'PASS' if passed else 'FAIL'}  (tolerance +-{int(TOL*100)}%)")
sys.exit(0 if passed else 1)
