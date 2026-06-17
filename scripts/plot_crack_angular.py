"""
plot_crack_angular.py
Compare simulated outgoing theta/phi distributions per crack against HFSS far-field theory.

Coordinate convention (from memory/project_coord_convention.md):
  normal_hat = x_world, theta_hat = y_world, phi_hat = z_world
  dir_out = sin(T)*cos(P)*x + sin(T)*sin(P)*y + cos(T)*z
  → HFSS_Theta = arccos(pz),  HFSS_Phi = atan2(py, px)
  Main transmission lobe: T=90°, P=0° → pure +x direction.

Usage:
    conda run -n bbrsim python scripts/plot_crack_angular.py [bbr_boundary_crossings.csv] [--iwt T] [--iwp P]

  --iwt / --iwp: fix the HFSS incoming wave angle (IWaveTheta / IWavePhi in degrees).
                 Use for gun-mode runs so theory matches the specific incidence angle.
                 Default: average over all (IWavePhi, IWaveTheta) combinations (Planck mode).
  Gun normal-incidence:  --iwt 180 --iwp 0
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("csv", nargs="?", default="build/output/bbr_boundary_crossings.csv")
parser.add_argument("--iwt", type=float, default=None,
                    help="Filter HFSS theory to this IWaveTheta [deg]. Default: average all.")
parser.add_argument("--iwp", type=float, default=None,
                    help="Filter HFSS theory to this IWavePhi [deg]. Default: average all.")
args = parser.parse_args()

CSV = args.csv
HFSS_BASE = "HFSSSimData"

CRACKS = {
    "crack1": "InfParallelPlate_crack1Rohan_500GHz",
    "crack2": "InfParallelPlate_crack2_500GHz",
}

# ── load simulated data ───────────────────────────────────────────────────────
df = pd.read_csv(CSV, on_bad_lines="skip", low_memory=False)
for col in ["px_post", "py_post", "pz_post", "n_reflect"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["px_post", "py_post", "pz_post"])

# Crack entry rows: photon arriving from world into vacuum_wg
# px_post > 0 → transmitted through crack (still heading in +x)
# px_post < 0 → reflected back (specular)
entering = df[(df["mat_pre"] != "vacuum_wg") & (df["mat_post"] == "vacuum_wg")]

def sim_angles(rows):
    """Convert transmitted crack-entry px/py/pz to HFSS (Theta, Phi) in degrees."""
    tx = rows[rows["px_post"] > 0].copy()
    pz = tx["pz_post"].values
    py = tx["py_post"].values
    px = tx["px_post"].values
    # clip for arccos numerical safety
    T = np.degrees(np.arccos(np.clip(pz, -1, 1)))   # 0=along z, 90=along x (forward)
    P = np.degrees(np.arctan2(py, px))              # azimuth around z-axis in xy-plane
    return T, P, len(tx)

# ── load HFSS far-field for one crack ────────────────────────────────────────
def load_hfss_farfield(dataset_id):
    """
    Load Ephi=0 and Ephi=1 far_field.csv for a crack dataset.
    Returns DataFrame with columns: IWavePhi, IWaveTheta, Phi, Theta, power
    where power = |rEphi|² + |rEtheta|²  (sum of Ephi=0 and Ephi=1 contributions).
    """
    rows = []
    for ephi in (0, 1):
        path = os.path.join(HFSS_BASE, f"{dataset_id}_Ephi={ephi}", "far_field.csv")
        d = pd.read_csv(path)
        d = d.rename(columns={"IWavePhi": "IWavePhi", "IWaveTheta": "IWaveTheta",
                               "Phi": "Phi", "Theta": "Theta"})
        d["power"] = (d["rEphi_real"]**2 + d["rEphi_imag"]**2 +
                      d["rEtheta_real"]**2 + d["rEtheta_imag"]**2)
        rows.append(d[["IWavePhi", "IWaveTheta", "Phi", "Theta", "power"]])
    combined = pd.concat(rows)
    # Sum Ephi=0 + Ephi=1 power at each (IWavePhi, IWaveTheta, Phi, Theta)
    return combined.groupby(["IWavePhi", "IWaveTheta", "Phi", "Theta"],
                             sort=False)["power"].sum().reset_index()

def theory_marginals(ff, iwt=None, iwp=None):
    """
    From a far-field DataFrame, compute theory marginal distributions for Theta and Phi.
    Weighting: sin(Theta) to match how BBRHFSSData weights for solid angle.

    iwt/iwp: if given, filter to the closest IWaveTheta/IWavePhi value before computing
             marginals (gun-mode: use specific incidence angle instead of averaging all).
    Returns: (theta_vals, theta_weights), (phi_vals, phi_weights)
    """
    ff = ff.copy()
    if iwt is not None:
        closest_t = ff["IWaveTheta"].unique()
        closest_t = closest_t[np.argmin(np.abs(closest_t - iwt))]
        ff = ff[ff["IWaveTheta"] == closest_t]
    if iwp is not None:
        closest_p = ff["IWavePhi"].unique()
        closest_p = closest_p[np.argmin(np.abs(closest_p - iwp))]
        ff = ff[ff["IWavePhi"] == closest_p]

    ff["sinT"] = np.sin(np.radians(ff["Theta"]))
    ff["w"] = ff["power"] * ff["sinT"]   # solid-angle-weighted power

    # Average over remaining incoming angle combinations
    n_incoming = ff.groupby(["IWavePhi", "IWaveTheta"]).ngroups
    ff_avg = ff.groupby(["Phi", "Theta"])["w"].sum().reset_index()
    ff_avg["w"] /= n_incoming

    # Marginal in Theta
    theta_m = ff_avg.groupby("Theta")["w"].sum().reset_index()
    theta_m = theta_m.sort_values("Theta")
    theta_w = theta_m["w"].values / theta_m["w"].sum()

    # Marginal in Phi
    phi_m = ff_avg.groupby("Phi")["w"].sum().reset_index()
    phi_m = phi_m.sort_values("Phi")
    phi_w = phi_m["w"].values / phi_m["w"].sum()

    return (theta_m["Theta"].values, theta_w), (phi_m["Phi"].values, phi_w)

# ── build plot ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
theory_label = (f"IWaveTheta={args.iwt}°, IWavePhi={args.iwp}°"
                if args.iwt is not None or args.iwp is not None
                else "averaged over all incoming angles")
fig.suptitle(f"Crack outgoing angle distributions: simulation vs HFSS theory\n"
             f"(main lobe: Theta=90°, Phi=0° → +x; theory: {theory_label})",
             fontsize=11)

# Offset by 2.5° so no HFSS Theta grid point (at 5° multiples) sits on a bin edge.
# Without offset, 6-dp rounding of pz causes arccos to return T slightly < the
# nominal value, mis-assigning HFSS-grid Theta values (10°, 50°, 100°, …) to
# the bin below.  Offset keeps every grid point ≥2.5° from the nearest edge.
T_BINS = np.arange(2.5, 183., 10.)   # edges: 2.5, 12.5, …, 172.5, 182.5
P_BINS = np.linspace(-90, 90, 19)    # 10° bins (Phi grid is not on 10° multiples)

for col, (label, dataset_id) in enumerate(CRACKS.items()):
    # Simulated
    crack_rows = entering[entering["vol_post"].str.contains(label.replace("crack", "crack"), na=False)]
    T_sim, P_sim, N_tx = sim_angles(crack_rows)
    N_ref = len(crack_rows) - N_tx
    print(f"{label}: {N_tx} transmitted, {N_ref} reflected")

    # HFSS theory
    ff = load_hfss_farfield(dataset_id)
    (theta_th, theta_w), (phi_th, phi_w) = theory_marginals(ff, iwt=args.iwt, iwp=args.iwp)

    # Theta panel
    ax = axes[0, col]
    if N_tx > 0:
        hn, he = np.histogram(T_sim, bins=T_BINS, density=False)
        hc = 0.5 * (he[:-1] + he[1:])
        ax.bar(hc, hn / hn.sum(), width=10, alpha=0.65, color="steelblue",
               label=f"Sim N={N_tx}")
    # Theory: bin discrete HFSS points using same convention as np.histogram ([lo, hi))
    th_bin_w, _ = np.histogram(theta_th, bins=T_BINS, weights=theta_w)
    th_bin_w = th_bin_w / th_bin_w.sum() if th_bin_w.sum() > 0 else th_bin_w
    tc = 0.5 * (T_BINS[:-1] + T_BINS[1:])
    ax.step(np.append(T_BINS[:-1], T_BINS[-1]),
            np.append(th_bin_w, th_bin_w[-1]),
            where="post", color="tomato", lw=2, label="HFSS theory")
    ax.axvline(90, color="gray", ls="--", lw=1, alpha=0.7, label="T=90° (forward)")
    ax.set_xlabel("HFSS Theta [deg]")
    ax.set_ylabel("Normalised counts")
    ax.set_title(f"{label} — Theta out\n(b={'50' if 'crack1' in label else '100'} µm)")
    ax.legend(fontsize=8)

    # Phi panel
    ax = axes[1, col]
    if N_tx > 0:
        hn, he = np.histogram(P_sim, bins=P_BINS, density=False)
        hc = 0.5 * (he[:-1] + he[1:])
        ax.bar(hc, hn / hn.sum(), width=10, alpha=0.65, color="darkorange",
               label=f"Sim N={N_tx}")
    ph_bin_w, _ = np.histogram(phi_th, bins=P_BINS, weights=phi_w)
    ph_bin_w = ph_bin_w / ph_bin_w.sum() if ph_bin_w.sum() > 0 else ph_bin_w
    ax.step(np.append(P_BINS[:-1], P_BINS[-1]),
            np.append(ph_bin_w, ph_bin_w[-1]),
            where="post", color="tomato", lw=2, label="HFSS theory")
    ax.axvline(0, color="gray", ls="--", lw=1, alpha=0.7, label="P=0° (forward)")
    ax.set_xlabel("HFSS Phi [deg]")
    ax.set_ylabel("Normalised counts")
    ax.set_title(f"{label} — Phi out")
    ax.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(os.path.dirname(CSV), "crack_angular_comparison.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
