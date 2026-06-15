"""
plot_test_output.py
Visualise test_output.csv: incoming Planck spectrum + outgoing photon distributions.

Usage:
    conda run -n bbrsim python scripts/plot_test_output.py [path/to/test_output.csv] [--temp T]
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── args ────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("csv", nargs="?", default="build/test_output.csv")
parser.add_argument("--temp", type=float, default=4.0, help="Emitter temperature [K]")
args = parser.parse_args()

df = pd.read_csv(args.csv, on_bad_lines="skip", low_memory=False)
for col in ["energy_eV", "n_reflect", "x_mm", "y_mm", "z_mm",
            "px_post", "py_post", "pz_post", "theta_in_deg", "phi_in_deg"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["energy_eV", "n_reflect"])

# Key events on (run_id, event_id) when the CSV carries run_id (multi-run
# sessions); fall back to event_id alone for older CSVs.
if "run_id" in df.columns:
    df["evt_key"] = list(zip(df["run_id"], df["event_id"]))
else:
    df["evt_key"] = df["event_id"]

T  = args.temp
kB = 8.617333e-5   # eV/K
kT = kB * T

print(f"Total rows      : {len(df)}")
print(f"Unique events   : {df['evt_key'].nunique()}")
print(f"Temperature     : {T} K  (kT = {kT*1e3:.3f} meV)")

# ── one row per photon for emission spectrum ─────────────────────────────────
first = df.groupby("evt_key", sort=False).first().reset_index()
E = first["energy_eV"].values

# ── Cu hits: photon arriving at a Cu surface (Drude names "Cu_RRR{N}_T{T}K") ──
cu_mask = df["mat_post"].str.startswith("Cu_RRR", na=False)
cu1 = df[cu_mask & (df["n_reflect"] == 1)]

# ── all Cu-surface hits (reflected or absorbed); the same row carries the
# incidence angle, and px/py/pz_post the outgoing direction if reflected ──────
cu_hits = df[cu_mask]

# ── photons entering crack volumes ───────────────────────────────────────────
cracks = df[df["mat_pre"] == "vacuum_wg"]

# ════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 10))
fig.suptitle(f"BBRsim test_output.csv  (T = {T} K,  N_events = {df['evt_key'].nunique():,})",
             fontsize=13, y=0.98)
gs = GridSpec(2, 3, figure=fig, hspace=0.40, wspace=0.35)

# ── Panel 1: Incoming energy spectrum vs Planck theory ────────────────────────
ax1 = fig.add_subplot(gs[0, :2])
u = E / kT
bins = np.linspace(0, 6, 31)
counts, edges = np.histogram(u, bins=bins)
centers = 0.5 * (edges[:-1] + edges[1:])
bin_w   = edges[1] - edges[0]
# Normalise histogram to integrate to 1
norm = counts.sum() * bin_w
ax1.bar(centers, counts / norm, width=bin_w, alpha=0.65,
        color="steelblue", label="Simulated")
# Planck photon-number spectrum: u²/(e^u-1), normalised
u_th  = np.linspace(0.01, 6, 400)
P_th  = u_th**2 / np.expm1(u_th)
norm_th = np.trapezoid(P_th, u_th)
ax1.plot(u_th, P_th / norm_th, "r-", lw=2, label=r"$u^2/(e^u-1)$ theory")
ax1.axvline(1.5936, color="gray", ls="--", lw=1, label="Peak u=1.594")
ax1.set_xlabel(r"$u = E / k_B T$")
ax1.set_ylabel("Probability density")
ax1.set_title("Incoming photon spectrum")
ax1.legend(fontsize=9)

# ── Panel 2: n_reflect histogram ──────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 2])
nr = df.groupby("evt_key")["n_reflect"].max().values
bins_n = np.arange(1, min(nr.max() + 2, 22))
cnts, edg = np.histogram(nr, bins=bins_n)
ax2.bar(0.5*(edg[:-1]+edg[1:]), cnts, width=1.0, alpha=0.75, color="darkorange")
ax2.set_yscale("log")
ax2.set_xlabel("n_reflect (reflections per track)")
ax2.set_ylabel("Counts")
ax2.set_title("Reflection count per photon")

# ── Panel 3: Cu first-hit position (y vs z) ──────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
if len(cu1) > 0:
    ax3.hexbin(cu1["z_mm"], cu1["y_mm"], gridsize=50, cmap="Blues",
               mincnt=1, linewidths=0.2)
    ax3.set_xlabel("z [mm]")
    ax3.set_ylabel("y [mm]")
    ax3.set_title(f"Cu first-hit position\n(N={len(cu1):,})")
else:
    ax3.text(0.5, 0.5, "No Cu first-hit events", ha="center", va="center",
             transform=ax3.transAxes)
    ax3.set_title("Cu first-hit position")

# ── Panel 4: Incidence angle (theta_in vs phi_in for Cu hits) ────────────────
ax4 = fig.add_subplot(gs[1, 1])
if len(cu_hits) > 0:
    ax4.hexbin(cu_hits["phi_in_deg"], cu_hits["theta_in_deg"],
               gridsize=40, cmap="Greens", mincnt=1, linewidths=0.2)
    ax4.set_xlabel("phi_in [deg]")
    ax4.set_ylabel("theta_in [deg]")
    ax4.set_title(f"Cu hits: incidence angle\n(N={len(cu_hits):,})")
else:
    ax4.text(0.5, 0.5, "No reflected events", ha="center", va="center",
             transform=ax4.transAxes)
    ax4.set_title("Cu hits: incidence angle")

# ── Panel 5: Crack entry positions ───────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2])
if len(cracks) > 0:
    c1e = cracks[cracks["vol_pre"].str.contains("crack1", na=False)]
    c2e = cracks[cracks["vol_pre"].str.contains("crack2", na=False)]
    if len(c1e): ax5.scatter(c1e["y_mm"], c1e["energy_eV"]*1e3, s=4,
                              alpha=0.5, label=f"crack1 (N={len(c1e)})", color="royalblue")
    if len(c2e): ax5.scatter(c2e["y_mm"], c2e["energy_eV"]*1e3, s=4,
                              alpha=0.5, label=f"crack2 (N={len(c2e)})", color="tomato")
    ax5.set_xlabel("y entry [mm]")
    ax5.set_ylabel("Energy [meV]")
    ax5.set_title("Crack entry photons")
    ax5.legend(fontsize=8, markerscale=3)
else:
    ax5.text(0.5, 0.5, "No crack events", ha="center", va="center",
             transform=ax5.transAxes)
    ax5.set_title("Crack entry photons")

import os
out = os.path.join(os.path.dirname(args.csv), "test_output_overview.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
