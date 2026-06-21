"""
check_nreflect.py
Validate the per-track reflection-count distribution from output/bbr.root.

PASS criteria:
  - n_reflect == 1 is the modal bin
  - counts are non-increasing for n in 1..10

Usage:
    conda run -n bbrsim python scripts/check_nreflect.py [path/to/bbr.root]
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"

df = load_crossings(PATH)
n = df["n_reflect"].values.astype(int)
bins = np.arange(1, min(n.max() + 2, 52))
counts, edges = np.histogram(n, bins=bins)

print(f"n_reflect range : [{n.min()}, {n.max()}]")
print(f"mean n_reflect  : {n.mean():.2f}")
print(f"modal bin       : {edges[np.argmax(counts)]:.0f}  (count={counts.max()})")

# In the single-wall test world max n_reflect == 1 is geometrically correct
# (a photon reflects once off Cu/crack, then exits; world-exit steps are
# fWorldBoundary and are not logged). "max > 1" is informational, not a gate.
ok_modal = int(edges[np.argmax(counts)]) == 1
window = counts[:min(10, len(counts))]
ok_mono = all(window[i] >= window[i + 1] for i in range(len(window) - 1))
passed = ok_modal and ok_mono

print(f"max n_reflect   : {n.max()}  (info only; 1 is expected here)")
print(f"mode = 1        : {'PASS' if ok_modal else 'FAIL'}")
print(f"monotone [1-10] : {'PASS' if ok_mono else 'FAIL'}")
print(f"RESULT          : {'PASS' if passed else 'FAIL'}")

fig, ax = plt.subplots(figsize=(8, 4))
centers = 0.5 * (edges[:-1] + edges[1:])
ax.bar(centers, counts, width=1.0, alpha=0.75, color="steelblue")
ax.set_yscale("log")
ax.set_xlabel("n_reflect (reflections per track)")
ax.set_ylabel("Counts (log scale)")
ax.set_title(f"Per-track reflection count  (N_crossings={len(df)})")
fig.tight_layout()
out = "nreflect_distribution.png"
fig.savefig(out, dpi=150)
print(f"Plot saved: {out}")

if not passed:
    sys.exit(1)
