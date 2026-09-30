"""
Validate that the BBRsim ROOT output (output/bbr.root) contains energies drawn from
the Planck photon-number spectrum at the given temperature.  The number spectrum
B ∝ ν²/(e^{hν/kT}−1) peaks at u = E/kT ≈ 1.5936.

Two rows, both required for PASS: the histogram peak (30 bins in u) within
[0.65, 1.35] of theory, and a Kolmogorov-Smirnov test of the first-crossing
energies against the photon-number CDF truncated to the emitter's fixed band
(4.14e-5 to 8.27e-2 eV), p > 0.01. The peak bin alone passes 10 K data analysed
at 9-12 K; the KS row rejects a 10 % temperature error.

Usage:
    conda run -n bbrsim python validation/check_planck_spectrum.py [path/to/bbr.root] [--temp T]
"""

import argparse
import sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bbrsim.io import load_crossings
from bbrsim import physics

parser = argparse.ArgumentParser()
parser.add_argument("csv", nargs="?", default="output/bbr.root")
parser.add_argument("--temp", type=float, default=4.0, help="Emitter temperature in K")
args = parser.parse_args()

CSV = args.csv
T   = args.temp

kT = physics.K_EV * T   # eV

df = load_crossings(CSV)
# n_reflect==1 is the first boundary crossing per track — captures emitted energy
data = df[df['n_reflect'] == 1]['energy_eV'].values
data = data[data > 0]
u    = data / kT

# Histogram in u
n_bins = 30
counts, edges = np.histogram(u, bins=n_bins, range=(0, 20))
centers = 0.5 * (edges[:-1] + edges[1:])

u_peak_obs = centers[np.argmax(counts)]
u_peak_theory = physics.PLANCK_PEAK_U

ratio = u_peak_obs / u_peak_theory
lo, hi = 0.65, 1.35
peak_ok = lo <= ratio <= hi

# KS against the photon-number CDF truncated to the emitter's band, built by a
# cumulative trapezoid on a fine log grid (the same p-values as a quad per sample).
E_MIN, E_MAX = 4.14e-5, 8.27e-2      # the emitter's fixed band [eV] (TestWorldPrimaryGeneratorAction)
grid = np.geomspace(E_MIN, E_MAX, 20001)
pdf = physics.planck_photon_number_pdf(grid, T)
cdf = np.concatenate([[0.], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
cdf /= cdf[-1]
ks_p = stats.kstest(data, lambda e: np.interp(e, grid, cdf)).pvalue if len(data) else 0.
ks_ok = ks_p > 0.01
passed = peak_ok and ks_ok

print(f"Events           : {len(data)}")
print(f"u_peak observed  : {u_peak_obs:.4f}")
print(f"u_peak theory    : {u_peak_theory:.4f}  (photon-number spectrum)")
print(f"ratio obs/theory : {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"KS p vs Planck   : {ks_p:.3g}  (expected > 0.01; CDF truncated to [{E_MIN}, {E_MAX}] eV)")
failed = [n for n, ok in (("peak", peak_ok), ("KS", ks_ok)) if not ok]
print(f"RESULT           : {'FAIL (' + ', '.join(failed) + ')' if failed else 'PASS'}")

# Plot
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(centers, counts, width=(edges[1]-edges[0]), alpha=0.7, label="Simulated")
ax.axvline(u_peak_obs, color="tab:blue", linestyle="--", label=f"Obs peak u={u_peak_obs:.3f}")
ax.axvline(u_peak_theory, color="tab:red", linestyle="-", label=f"Theory peak u={u_peak_theory:.4f}")
ax.set_xlabel("u = E / kT")
ax.set_ylabel("Counts")
ax.set_title(f"Planck photon-number spectrum at T={T:.4g} K  (ratio={ratio:.3f})")
ax.legend()
fig.tight_layout()
out = "planck_spectrum_check.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
