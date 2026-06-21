"""
check_crack_ratio.py
Compare crack2/crack1 event-rate ratio to the expected geometric aperture ratio.
Both cracks have the same HFSS transmittance (~52%), so the rate ratio should
equal the opening-area ratio. PASS if within 3 sigma_Poisson of expected.

Usage:
    conda run -n bbrsim python scripts/check_crack_ratio.py [path/to/bbr.root]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings
from bbrsim import select

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"

# Aperture areas [mm^2]
A1 = 2 * 5.1 * 2 * 0.026   # crack1
A2 = 2 * 5.1 * 2 * 0.051   # crack2
expected_ratio = A2 / A1

df = load_crossings(PATH)
wg = select.crack_crossings(df)
c1 = wg[wg["vol_pre"].str.contains("crack1", na=False)]
c2 = wg[wg["vol_pre"].str.contains("crack2", na=False)]

N1, N2 = len(c1), len(c2)
print(f"crack1 events : {N1}")
print(f"crack2 events : {N2}")
print(f"aperture A1   : {A1:.4f} mm^2")
print(f"aperture A2   : {A2:.4f} mm^2")
print(f"expected ratio (A2/A1) : {expected_ratio:.3f}")

if N1 == 0:
    print("ERROR: no crack1 events - cannot compute ratio")
    sys.exit(1)

obs_ratio = N2 / N1
sigma_ratio = obs_ratio * np.sqrt(1 / N2 + 1 / N1) if N2 > 0 else float("inf")
n_sigma = abs(obs_ratio - expected_ratio) / sigma_ratio
passed = n_sigma < 3.0

print(f"observed ratio N2/N1   : {obs_ratio:.3f} +/- {sigma_ratio:.3f}")
print(f"deviation (sigma)      : {n_sigma:.2f}  (threshold < 3)")
print(f"RESULT                 : {'PASS' if passed else 'FAIL'}")

if not passed:
    sys.exit(1)
