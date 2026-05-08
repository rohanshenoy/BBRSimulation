"""
Validate that planck_output.csv contains energies drawn from the Planck photon-number
spectrum at 4 K.  The number spectrum B ∝ ν²/(e^{hν/kT}−1) peaks at u = E/kT ≈ 1.5936.

Usage:
    conda run -n bbrsim python scripts/check_planck_spectrum.py [path/to/planck_output.csv]
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "build/planck_output.csv"

# Physical constants
k_eV = 8.6173e-5  # eV/K
T    = 4.0         # K
kT   = k_eV * T   # eV

data = np.loadtxt(CSV, skiprows=1)
u    = data / kT

# Histogram in u
n_bins = 100
counts, edges = np.histogram(u, bins=n_bins, range=(0, 20))
centers = 0.5 * (edges[:-1] + edges[1:])

u_peak_obs = centers[np.argmax(counts)]
u_peak_theory = 1.5936

ratio = u_peak_obs / u_peak_theory
lo, hi = 0.65, 1.35
passed = lo <= ratio <= hi

print(f"Events           : {len(data)}")
print(f"u_peak observed  : {u_peak_obs:.4f}")
print(f"u_peak theory    : {u_peak_theory:.4f}  (photon-number spectrum)")
print(f"ratio obs/theory : {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT           : {'PASS' if passed else 'FAIL'}")

# Plot
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(centers, counts, width=(edges[1]-edges[0]), alpha=0.7, label="Simulated")
ax.axvline(u_peak_obs, color="tab:blue", linestyle="--", label=f"Obs peak u={u_peak_obs:.3f}")
ax.axvline(u_peak_theory, color="tab:red", linestyle="-", label=f"Theory peak u={u_peak_theory:.4f}")
ax.set_xlabel("u = E / kT")
ax.set_ylabel("Counts")
ax.set_title(f"Planck photon-number spectrum at T=4 K  (ratio={ratio:.3f})")
ax.legend()
fig.tight_layout()
out = "planck_spectrum_check.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
