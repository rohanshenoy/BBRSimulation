"""
check_crack_wall_reflection.py
Validate the reflection normal when a photon inside a vacuum_wg crack strikes
the crack's side wall (crack daughter -> Cu slab mother). Run on the output of
crack_wall.mac (gun inside crack1 heading mostly +z, so every photon hits the
z-wall at z = +0.026 mm).

Expected physics: specular reflection about the WALL normal (z):
    pz_post = -pz_pre,   px_post = px_pre,   py_post = py_pre.
The break this catches: taking the normal from the entered (slab) solid,
which for an interior point returns the slab's x-face normal, so px flips
and pz does not — the photon then continues into the copper.

PASS if every crack->Cu BBRReflect row flips pz and preserves px, py
(|delta| < 1e-9), and at least one such row exists.

Usage:
    conda run -n bbrsim python scripts/check_crack_wall_reflection.py [path/to/bbr.root]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim.io import load_crossings

PATH = sys.argv[1] if len(sys.argv) > 1 else "build/output/bbr.root"
TOL = 1e-9

df = load_crossings(PATH)
wall = df[(df["mat_pre"] == "vacuum_wg")
          & df["mat_post"].str.startswith("Cu_RRR", na=False)
          & (df["status"] == "BBRReflect")]

n = len(wall)
print(f"crack -> Cu wall reflections : {n}")
if n == 0:
    print("no crack-wall reflections found — did you run crack_wall.mac?")
    print("RESULT: FAIL")
    sys.exit(1)

pz_flipped  = np.abs(wall["pz_post"].values + wall["pz_pre"].values) < TOL
px_kept     = np.abs(wall["px_post"].values - wall["px_pre"].values) < TOL
py_kept     = np.abs(wall["py_post"].values - wall["py_pre"].values) < TOL
px_flipped  = np.abs(wall["px_post"].values + wall["px_pre"].values) < TOL

print(f"pz flipped (correct, z-wall) : {pz_flipped.sum()}")
print(f"px preserved                 : {px_kept.sum()}")
print(f"py preserved                 : {py_kept.sum()}")
print(f"px flipped (WRONG, slab face): {px_flipped.sum()}")

passed = bool(pz_flipped.all() and px_kept.all() and py_kept.all())
print(f"RESULT: {'PASS' if passed else 'FAIL'}")
sys.exit(0 if passed else 1)
