"""
Validate that test_output.csv contains energies drawn from the Planck photon-number
spectrum at 4 K.  The number spectrum B ∝ ν²/(e^{hν/kT}−1) peaks at u = E/kT ≈ 1.5936.

Usage:
    conda run -n bbrsim python scripts/check_planck_spectrum.py [path/to/test_output.csv]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "build/test_output.csv"

# Physical constants
k_eV = 8.6173e-5  # eV/K
T    = 4.0         # K
kT   = k_eV * T   # eV

df   = pd.read_csv(CSV, on_bad_lines='skip')
df['energy_eV'] = pd.to_numeric(df['energy_eV'], errors='coerce')
df['n_reflect'] = pd.to_numeric(df['n_reflect'], errors='coerce')
df = df.dropna(subset=['energy_eV', 'n_reflect'])
# n_reflect==1 selects the first boundary event per photon track = emitted energy.
# If n_reflect is a global counter (not per-track), fall back to first row per event_id.
if (df['n_reflect'] == 1).sum() < 2:
    first = df.groupby('event_id', sort=False).first().reset_index()
    data = first['energy_eV'].dropna().values
else:
    data = df[df['n_reflect'] == 1]['energy_eV'].values
# Drop unphysical energies (corrupted rows)
data = data[data > 0]
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
