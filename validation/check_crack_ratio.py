"""
check_crack_ratio.py
Compare crack2/crack1 event-rate ratio to the expected geometric aperture ratio.
The check counts crack ENTRIES (independent of the HFSS transmittance), so the
rate ratio should equal the opening-area ratio when both cracks see the same
illumination. PASS if within 3 sigma_Poisson of expected.

Input: the second run of validation/G4Macros/Validation_CrackTransmit.mac,
output/bbr_ratio.root (the regression runner reads it from
<build>/regression/transmit/output/): 200 000 Planck photons at 4 K from a thin
slab just off the Cu face, centred between the two cracks so both see the same
illumination, which gives ~640 / ~1290 entries. Any Planck run that
illuminates both cracks alike works too, e.g. planck_10K.mac (1M events,
~70 / ~180); planck.mac gives too few, and a fixed-gun run at one crack has no
ratio to test.

Usage (the path defaults to output/bbr.root):
    conda run -n bbrsim python <repo>/validation/check_crack_ratio.py output/bbr_ratio.root
"""

import sys

import numpy as np

from bbrsim.io import load_crossings
from bbrsim import select

PATH = sys.argv[1] if len(sys.argv) > 1 else "output/bbr.root"

# Aperture areas [mm^2]
A1 = 2 * 5.1 * 2 * 0.026   # crack1
A2 = 2 * 5.1 * 2 * 0.051   # crack2
expected_ratio = A2 / A1

df = load_crossings(PATH)
wg = select.crack_crossings(df)
c1 = wg[wg["vol_post"].str.contains("crack1", na=False)]
c2 = wg[wg["vol_post"].str.contains("crack2", na=False)]

N1, N2 = len(c1), len(c2)
print(f"crack1 events : {N1}")
print(f"crack2 events : {N2}")
print(f"aperture A1   : {A1:.4f} mm^2")
print(f"aperture A2   : {A2:.4f} mm^2")
print(f"expected ratio (A2/A1) : {expected_ratio:.3f}")

if N1 == 0 or N2 == 0:
    print(f"ERROR: no {'crack1' if N1 == 0 else 'crack2'} events - cannot compute "
          "ratio (needs a Planck run that illuminates both cracks: the second run of "
          "Validation_CrackTransmit.mac, output/bbr_ratio.root)")
    sys.exit(1)

obs_ratio = N2 / N1
sigma_ratio = obs_ratio * np.sqrt(1 / N2 + 1 / N1)
n_sigma = abs(obs_ratio - expected_ratio) / sigma_ratio
passed = n_sigma < 3.0

print(f"observed ratio N2/N1   : {obs_ratio:.3f} +/- {sigma_ratio:.3f}")
print(f"deviation (sigma)      : {n_sigma:.2f}  (threshold < 3)")
print(f"RESULT                 : {'PASS' if passed else 'FAIL'}")

if not passed:
    sys.exit(1)
