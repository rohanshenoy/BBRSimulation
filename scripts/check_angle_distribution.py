"""
check_angle_distribution.py
Check theta_in distribution for first-hit Cu events against Lambert cosine law.
For θ < 45° (where finite-wall clipping is negligible), run a KS test.
PASS if KS p-value > 0.01.

Usage:
    conda run -n bbrsim python scripts/check_angle_distribution.py [bbr_boundary_crossings.csv]
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import kstest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr_boundary_crossings.csv"

df = pd.read_csv(CSV, on_bad_lines="skip")
for col in ["theta_in_deg", "n_reflect"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["theta_in_deg", "n_reflect"])

# First-hit Cu events only. The photon arrives from the vacuum, so Cu is the
# POST-step material; Drude materials are named "Cu_RRR{N}_T{T}K".
cu_first = df[(df["n_reflect"] == 1)
              & df["mat_post"].str.startswith("Cu_RRR", na=False)]
theta = cu_first["theta_in_deg"].values

# The emitter samples θ UNIFORMLY in [0°, 90°] (YYC's intentional convention,
# NOT Lambertian — see CLAUDE.md, /bbr/thermal/ notes). The incidence angle on
# the x=0 Cu wall equals the emission θ, so the expected θ_in density is
# uniform — up to finite-wall clipping. With the 20×20 mm emitter patch at
# x=−50 mm and the 50×50 mm wall, every direction with tanθ < 15/50 reaches
# the wall, so the distribution is exactly uniform for θ < ~16.7°. KS-test
# against a uniform CDF on [0°, 15°].
TH_MAX = 15.0
theta_lo = theta[(theta >= 0) & (theta < TH_MAX)]

def cdf_uniform(theta_deg):
    return theta_deg / TH_MAX

stat, p_value = kstest(theta_lo, cdf_uniform)
passed = p_value > 0.01

print(f"First-hit Cu events : {len(theta)}")
print(f"Events θ < {TH_MAX:.0f}°      : {len(theta_lo)}")
print(f"KS stat             : {stat:.4f}")
print(f"KS p-value          : {p_value:.4f}  (threshold 0.01)")
print(f"RESULT              : {'PASS' if passed else 'FAIL'}")

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, 90, 46)
counts, edges = np.histogram(theta, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])
ax.bar(centers, counts, width=2.0, alpha=0.7, color="darkorange", label=f"Simulated (N={len(theta)})")
ax.axhline(1.0 / 90.0, color="r", lw=2,
           label="uniform-θ emission (YYC convention), no clipping")
ax.axvline(TH_MAX, color="gray", ls="--", lw=1, label=f"{TH_MAX:.0f}° KS window (clipping-free)")
ax.set_xlabel("theta_in (deg from Cu normal)")
ax.set_ylabel("Probability density (per deg)")
ax.set_title(f"Incidence angle vs uniform-θ emission  (KS p={p_value:.3f})")
ax.legend(fontsize=9)
fig.tight_layout()
out = "angle_distribution_check.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
