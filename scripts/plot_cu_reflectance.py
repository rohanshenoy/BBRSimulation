"""
plot_cu_reflectance.py
Plot copper reflectance (R) and absorptance (D = 1-R) vs frequency.
Three panels:
  1. Reflectance R — full Drude for OFHC/OF/HP Cu + tabulated data
  2. Absorptance D = 1-R — same data
  3. Temperature dependence — full Drude for OFHC Cu (RRR=100) at 4, 20, 77, 300 K

All Cu curves use the full Drude model parameterized by RRR (Griffiths §9.4).
Hagen-Rubens is shown only for OFHC_Cu as a low-frequency reference.
σ_phonon ~ 1/T is the simple power-law approximation valid for T > ~50 K;
below ~50 K umklapp scattering dominates — the simple formula is not used
(σ_DC = RRR × σ_RT is used instead). This only matters for warm shield layers.

Data sources:
  - Full Drude model for Cu at RRR=100 (OFHC_Cu), RRR=3 (OF_Cu), RRR=6 (HP_Cu)
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
m_e  = 9.109e-31          # kg
n_e  = 8.49e28            # free electrons/m³ (copper)
e_C  = 1.602e-19          # C

# ── Hagen-Rubens model (low-frequency limit only) ────────────────────────────
def hagen_rubens_R(freq_Hz, sigma_SI):
    omega = 2.0 * np.pi * freq_Hz
    D = 2.0 * np.sqrt(2.0 * eps0 * omega / sigma_SI)
    return 1.0 - np.clip(D, 0, 1)

# ── Full Drude model (Griffiths §9.4) ────────────────────────────────────────
sigma_RT = 5.96e7   # S/m — universal for all Cu grades at 273 K
# RRR is the primary user parameter. σ_DC = RRR × σ_RT at 4K (impurity dominated).
# σ_phonon ~ 1/T (simple power law) is only valid above ~50K.
# Below ~50K umklapp scattering makes the phonon term more complex (Bloch-Grüneisen)
# but it is negligible vs the impurity term at 4K anyway.

def sigma_drude(T_K, RRR):
    """DC conductivity via Matthiessen's rule.
    Linear phonon approximation valid for T >= 50 K.
    Below 50 K phonons are frozen; sigma_DC = RRR * sigma_RT.
    """
    sigma_imp = RRR * sigma_RT
    if T_K >= 50.0:
        sigma_ph = sigma_RT * 273.0 / T_K
        return 1.0 / (1.0/sigma_imp + 1.0/sigma_ph)
    return sigma_imp

def drude_R(freq_Hz, sigma_DC):
    """Normal-incidence reflectance from the full complex Drude model.
    sigma_DC: DC conductivity [S/m]; tau derived from Drude formula.

    σ(ω) = σ_DC/(1−iωτ) is kept complex: ε̃ = 1 + iσ/(ε₀ω), ñ = √ε̃,
    R = |(ñ−1)/(ñ+1)|². Im σ supplies the plasma term in Re ε̃, which keeps
    R near 1 in the relaxation regime (D ≈ 2/(ωp·τ), flat in frequency).
    """
    tau   = sigma_DC * m_e / (n_e * e_C**2)
    omega = 2.0 * np.pi * freq_Hz
    sigma_ac = sigma_DC / (1.0 - 1j * omega * tau)
    eps_t = 1.0 + 1j * sigma_ac / (eps0 * omega)
    n_t   = np.sqrt(eps_t)
    R     = np.abs((n_t - 1.0) / (n_t + 1.0)) ** 2
    return np.clip(R, 0.0, 1.0)

# Drude materials — RRR is the primary parameter; σ_DC = RRR × σ_RT at 4K.
# Named aliases match GetCopperByName() in BBRMaterials.hh.
materials_drude = {
    "OFHC_Cu (RRR=100, Drude, 4K)": {"rrr": 100, "color": "steelblue",   "ls": "-"},
    "OF_Cu   (RRR=3,   Drude, 4K)": {"rrr":   3, "color": "darkorange",  "ls": "--"},
    "HP_Cu   (RRR=6,   Drude, 4K)": {"rrr":   6, "color": "forestgreen", "ls": "-."},
}

# ── load tabulated data ───────────────────────────────────────────────────────
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "materials")

g4ir   = pd.read_csv(os.path.join(data_dir, "cu_g4ir_reflectivity.csv"))
palik  = pd.read_csv(os.path.join(data_dir, "cu_palik_optical_constants.csv"))
serov  = pd.read_csv(os.path.join(data_dir, "cu_serov_reference_points.csv"))

# ── frequency grid for H-R curves ────────────────────────────────────────────
freq_GHz = np.logspace(-1, 5, 600)   # 0.1 GHz to 100 THz
freq_Hz  = freq_GHz * 1e9

# ── figure: three panels ─────────────────────────────────────────────────────
fig, (ax_R, ax_D, ax_T) = plt.subplots(1, 3, figsize=(20, 6))
fig.suptitle("Copper reflectance & absorptance vs frequency", fontsize=13)

for ax, qty, ylabel, ylim, yscale in [
    (ax_R, "R", "Reflectance  R",       (0.93, 1.002), "log"),
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

    # Drude curves — one per RRR alias; this is what BBRMaterials.hh computes
    for label, cfg in materials_drude.items():
        sig  = sigma_drude(4.0, cfg["rrr"])
        vals = drude_R(freq_Hz, sig)
        y    = vals if qty == "R" else 1 - vals
        ax.plot(freq_GHz, y, color=cfg["color"], ls=cfg["ls"], lw=2, label=label)

    # Hagen-Rubens for OFHC_Cu — shown as dashed reference to mark where it diverges
    sigma_ofhc = 100 * sigma_RT
    hr_vals = hagen_rubens_R(freq_Hz, sigma_ofhc)
    y_hr = hr_vals if qty == "R" else 1 - hr_vals
    ax.plot(freq_GHz, y_hr, color="steelblue", ls=":", lw=1.5, alpha=0.5,
            label="OFHC_Cu Hagen-Rubens (low-f limit, ref only)")

    # G4 REFLECTIVITY table: the 24 log-spaced points stored by BuildDrudeMaterial
    # (same Drude formula evaluated at discrete energies — what Geant4 actually
    # interpolates). 10 GHz–20 THz, matching the Planck-emitter CDF range.
    N_tab      = 24
    E_tab_eV   = np.exp(np.linspace(np.log(4.14e-5), np.log(8.27e-2), N_tab))  # eV
    h_eVs_loc  = 4.13566769692e-15
    nu_tab     = E_tab_eV / h_eVs_loc
    fGHz_tab   = nu_tab / 1e9
    R_tab      = drude_R(nu_tab * 1e9, sigma_ofhc)
    y_tab = R_tab if qty == "R" else 1 - R_tab
    ax.plot(fGHz_tab, y_tab, "P", color="steelblue", ms=6, alpha=0.9, zorder=6,
            markeredgecolor="navy", markeredgewidth=0.6,
            label="G4 Drude table (24 pts, Cu_RRR100_T4K in BBRMaterials)")

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

# ── panel 3: temperature dependence (Drude, OFHC Cu RRR=100) ─────────────────
temps = [4, 20, 77, 300]
cmap  = plt.cm.plasma
tcolors = [cmap(0.05), cmap(0.3), cmap(0.6), cmap(0.9)]

ax_T.set_xscale("log")
ax_T.set_yscale("log")
ax_T.set_xlabel("Frequency [GHz]")
ax_T.set_ylabel("Absorptance  D = 1−R")
ax_T.set_xlim(0.1, 1e5)
ax_T.set_ylim(1e-5, 0.1)
ax_T.grid(True, which="both", ls=":", alpha=0.4)
ax_T.axvspan(50, 20000, alpha=0.06, color="gray", label="BBRsim range")
ax_T.set_title("OFHC Cu (RRR=100) — temperature dependence", fontsize=10)

for T, col in zip(temps, tcolors):
    sig = sigma_drude(T, RRR=100)
    D_drude = 1.0 - drude_R(freq_Hz, sig)
    tau_ps  = sig * m_e / (n_e * e_C**2) * 1e12
    f_break = 1.0 / (2.0 * np.pi * sig * m_e / (n_e * e_C**2)) / 1e9
    linestyle = "-" if T >= 50 else "--"  # solid above ~50K (1/T valid), dashed below
    ax_T.plot(freq_GHz, D_drude, color=col, lw=2, ls=linestyle,
              label=f"T = {T} K  (σ={sig:.2e} S/m, f_break={f_break:.0f} GHz)")

# Mark the ~50K boundary below which σ_phonon ~ 1/T breaks down (umklapp regime)
ax_T.axhline(y=0, alpha=0)  # invisible anchor
ax_T.annotate(
    "Below ~50 K: σ_phonon ~ 1/T\nbreaks down (umklapp/Bloch-Grüneisen).\n"
    "σ_DC = RRR × σ_RT used instead.",
    xy=(0.3, 0.03), xycoords="axes fraction",
    fontsize=7, color="gray",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7),
)

# Serov 4K reference points on temperature panel
serov_colors = {"OF_Cu": "darkorange", "HP_Cu": "forestgreen"}
for _, row in serov.iterrows():
    ax_T.scatter(row["freq_GHz"], row["D"], s=80, zorder=5,
                 color=serov_colors.get(row["material"], "red"),
                 edgecolors="black", linewidths=0.8,
                 label=f"{row['material']} @ {row['freq_GHz']:.0f} GHz (Serov 4K)")

ax_T.legend(fontsize=7.5, loc="upper left")

# ── frequency tick labels ─────────────────────────────────────────────────────
for ax in (ax_R, ax_D, ax_T):
    ax.set_xticks([1, 10, 100, 1000, 10000, 100000])
    ax.set_xticklabels(["1 GHz", "10", "100", "1 THz", "10", "100 THz"])

plt.tight_layout()
os.makedirs(os.path.dirname(args.out), exist_ok=True)
fig.savefig(args.out, dpi=150, bbox_inches="tight")
print(f"Saved: {args.out}")
