"""
check_angle_distribution.py
Check theta_in distribution for first-hit Cu events against Lambert cosine law.
For θ < 45° (where finite-wall clipping is negligible), run a KS test.
PASS if KS p-value > 0.01.

Usage:
    conda run -n bbrsim python scripts/check_angle_distribution.py [test_output.csv]
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import kstest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "build/test_output.csv"

df = pd.read_csv(CSV, on_bad_lines="skip")
for col in ["theta_in_deg", "n_reflect"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["theta_in_deg", "n_reflect"])

# First-hit Cu events only
cu_first = df[(df["n_reflect"] == 1) & (df["mat_pre"] == "OFHC_Cu")]
theta = cu_first["theta_in_deg"].values

# KS test against sin(2θ) CDF for θ ∈ [0, 45°] only (avoid wall-clipping bias)
theta_lo = theta[(theta >= 0) & (theta < 45)]

# CDF of sin(2θ) on [0, 45°]: F(θ) = (1 - cos(2θ)) / (1 - cos(90°)) = (1 - cos(2θ)) / 1
def cdf_lambert(theta_deg):
    t = np.radians(theta_deg)
    return (1 - np.cos(2 * t)) / 1.0  # normalised to [0,45°] window

stat, p_value = kstest(theta_lo, cdf_lambert)
passed = p_value > 0.01

print(f"First-hit Cu events : {len(theta)}")
print(f"Events θ < 45°      : {len(theta_lo)}")
print(f"KS stat             : {stat:.4f}")
print(f"KS p-value          : {p_value:.4f}  (threshold 0.01)")
print(f"RESULT              : {'PASS' if passed else 'FAIL'}")

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, 90, 46)
counts, edges = np.histogram(theta, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])
ax.bar(centers, counts, width=2.0, alpha=0.7, color="darkorange", label=f"Simulated (N={len(theta)})")
t_th = np.linspace(0, 90, 300)
pdf_th = np.sin(np.radians(2 * t_th)) * np.pi / 180  # convert to per-deg
# Normalise over [0,90]: integral = 1 by construction
ax.plot(t_th, pdf_th, "r-", lw=2, label="sin(2θ) Lambert")
ax.axvline(45, color="gray", ls="--", lw=1, label="45° KS limit")
ax.set_xlabel("theta_in (deg from Cu normal)")
ax.set_ylabel("Probability density (per deg)")
ax.set_title(f"Incidence angle vs Lambert cosine law  (KS p={p_value:.3f})")
ax.legend(fontsize=9)
fig.tight_layout()
out = "angle_distribution_check.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
