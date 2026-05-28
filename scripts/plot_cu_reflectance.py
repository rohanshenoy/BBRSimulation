"""
plot_cu_reflectance.py
Plot copper reflectance (R) and absorptance (D = 1-R) vs frequency,
overlaying all available data sources:
  - Hagen-Rubens analytical curves (OFHC_Cu, OF_Cu, HP_Cu)
  - Geant4 IR reflectivity table (yyc / Geant4_copper_IR_reflectivity.ods)
  - Palik Handbook of Optical Constants, Vol. 1 Table 1 (room temperature)
  - Serov et al. (2016) cryogenic reference points at 4 K

Usage:
    conda run -n bbrsim python scripts/plot_cu_reflectance.py [--out path/to/out.png]
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--out", default="build/cu_reflectance_plots.png")
args = parser.parse_args()

# ── physical constants ────────────────────────────────────────────────────────
eps0 = 8.8541878128e-12   # F/m

# ── Hagen-Rubens model ────────────────────────────────────────────────────────
def hagen_rubens_R(freq_Hz, sigma_SI):
    omega = 2.0 * np.pi * freq_Hz
    D = 2.0 * np.sqrt(2.0 * eps0 * omega / sigma_SI)
    return 1.0 - np.clip(D, 0, 1)

# Conductivities from BBRMaterials.hh
materials_HR = {
    "OFHC_Cu (RRR=100, Hagen-Rubens)": {"sigma": 5.96e9,  "color": "steelblue",   "ls": "-"},
    "OF_Cu (Serov 2016, Hagen-Rubens)": {"sigma": 1.786e8, "color": "darkorange",  "ls": "--"},
    "HP_Cu (Serov 2016, Hagen-Rubens)": {"sigma": 3.383e8, "color": "forestgreen", "ls": "-."},
}

# ── load tabulated data ───────────────────────────────────────────────────────
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

g4ir   = pd.read_csv(os.path.join(data_dir, "cu_g4ir_reflectivity.csv"))
palik  = pd.read_csv(os.path.join(data_dir, "cu_palik_optical_constants.csv"))
serov  = pd.read_csv(os.path.join(data_dir, "cu_serov_reference_points.csv"))

# ── frequency grid for H-R curves ────────────────────────────────────────────
freq_GHz = np.logspace(-1, 5, 600)   # 0.1 GHz to 100 THz
freq_Hz  = freq_GHz * 1e9

# ── figure: two panels (R and D = 1-R) ───────────────────────────────────────
fig, (ax_R, ax_D) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Copper reflectance & absorptance vs frequency", fontsize=13)

for ax, qty, ylabel, ylim, yscale in [
    (ax_R, "R", "Reflectance  R",       (0.94, 1.002), "linear"),
    (ax_D, "D", "Absorptance  D = 1−R", (1e-5, 0.1),   "log"),
]:
    ax.set_xscale("log")
    ax.set_yscale(yscale)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.1, 1e5)
    ax.set_ylim(*ylim)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.axvspan(50, 20000, alpha=0.06, color="gray", label="BBRsim range (50 GHz–20 THz)")

    # Hagen-Rubens curves
    for label, cfg in materials_HR.items():
        vals = hagen_rubens_R(freq_Hz, cfg["sigma"])
        y = vals if qty == "R" else 1 - vals
        ax.plot(freq_GHz, y, color=cfg["color"], ls=cfg["ls"], lw=2, label=label)

    # Geant4 IR reflectivity table
    y_g4 = g4ir["R"] if qty == "R" else g4ir["D"]
    ax.plot(g4ir["freq_GHz"], y_g4, "k^", ms=4, alpha=0.7,
            label="G4 IR table (Geant4_copper_IR_reflectivity.ods)")

    # Palik handbook (room temperature) — restrict to free-electron regime (k > 1).
    # The 1-2 THz band in Palik Vol.1 has anomalously low k (~0.1), likely a
    # measurement-gap artefact; those points are excluded.
    palik_plot = palik[(palik["freq_GHz"] >= 0.1) & (palik["freq_GHz"] <= 1e5)
                       & (palik["k"] > 1.0)]
    y_p = palik_plot["R"] if qty == "R" else palik_plot["D"]
    ax.plot(palik_plot["freq_GHz"], y_p, "ms", ms=4, alpha=0.6,
            label="Palik Handbook Vol.1, Table 1 (room temp, k>1)")

    # Serov cryogenic reference points
    serov_colors = {"OF_Cu": "darkorange", "HP_Cu": "forestgreen"}
    for _, row in serov.iterrows():
        y_s = row["R"] if qty == "R" else row["D"]
        ax.scatter(row["freq_GHz"], y_s, s=80, zorder=5,
                   color=serov_colors.get(row["material"], "red"),
                   edgecolors="black", linewidths=0.8,
                   label=f"{row['material']} @ {row['freq_GHz']:.0f} GHz ({row['source']})")

    ax.legend(fontsize=7.5, loc="lower left" if qty == "R" else "upper left")

# ── frequency tick labels ─────────────────────────────────────────────────────
for ax in (ax_R, ax_D):
    ax.set_xticks([1, 10, 100, 1000, 10000, 100000])
    ax.set_xticklabels(["1 GHz", "10", "100", "1 THz", "10", "100 THz"])

plt.tight_layout()
os.makedirs(os.path.dirname(args.out), exist_ok=True)
fig.savefig(args.out, dpi=150, bbox_inches="tight")
print(f"Saved: {args.out}")
