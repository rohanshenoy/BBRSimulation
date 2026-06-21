"""
check_angle_distribution.py
Check theta_in distribution for first-hit Cu events against the emitter's
uniform-in-theta convention. KS-test the clipping-free window (theta < 15 deg).
PASS if KS p-value > 0.01.

Usage:
    conda run -n bbrsim python scripts/check_angle_distribution.py [path/to/bbr.root]
"""

import os
import sys

import numpy as np
from scipy.stats import kstest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings
from bbrsim import select

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"

df = load_crossings(PATH)
theta = select.first_hit_cu(df)["theta_in_deg"].values

# The emitter samples theta UNIFORMLY in [0, 90] deg (YYC's intentional
# convention, NOT Lambertian). The incidence angle on the x=0 Cu wall equals
# the emission theta, so theta_in is uniform up to finite-wall clipping. With
# the 20x20 mm patch at x=-50 mm and the 50x50 mm wall, every direction with
# tan(theta) < 15/50 reaches the wall, so the distribution is uniform for
# theta < ~16.7 deg. KS-test against a uniform CDF on [0, 15] deg.
TH_MAX = 15.0
theta_lo = theta[(theta >= 0) & (theta < TH_MAX)]


def cdf_uniform(theta_deg):
    return theta_deg / TH_MAX


stat, p_value = kstest(theta_lo, cdf_uniform)
passed = p_value > 0.01

print(f"First-hit Cu events : {len(theta)}")
print(f"Events theta < {TH_MAX:.0f}   : {len(theta_lo)}")
print(f"KS stat             : {stat:.4f}")
print(f"KS p-value          : {p_value:.4f}  (threshold 0.01)")
print(f"RESULT              : {'PASS' if passed else 'FAIL'}")

fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, 90, 46)
counts, edges = np.histogram(theta, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])
ax.bar(centers, counts, width=2.0, alpha=0.7, color="darkorange",
       label=f"Simulated (N={len(theta)})")
ax.axhline(1.0 / 90.0, color="r", lw=2,
           label="uniform-theta emission (YYC convention), no clipping")
ax.axvline(TH_MAX, color="gray", ls="--", lw=1, label=f"{TH_MAX:.0f} deg KS window")
ax.set_xlabel("theta_in (deg from Cu normal)")
ax.set_ylabel("Probability density (per deg)")
ax.set_title(f"Incidence angle vs uniform-theta emission  (KS p={p_value:.3f})")
ax.legend(fontsize=9)
fig.tight_layout()
out = "angle_distribution_check.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
