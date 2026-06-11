"""
plot_brief_figures.py
Generate figures for copper_reflectance_brief.tex.

Outputs (same directory as this script):
  fig_drude_scaling.pdf   -- universal Drude D/D_HR curve + AC conductivity
  fig_temp_dependence.pdf -- D vs frequency for T = 4, 20, 77, 300 K

Usage:
    conda run -n bbrsim python docs/physics/plot_brief_figures.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.family':      'serif',
    'font.size':        11,
    'axes.labelsize':   11,
    'axes.titlesize':   11,
    'legend.fontsize':   9,
    'xtick.labelsize':   9,
    'ytick.labelsize':   9,
    'figure.dpi':       150,
})

# ── physical constants ────────────────────────────────────────────────────────
eps0     = 8.8541878128e-12   # F/m
m_e      = 9.109e-31          # kg
n_e      = 8.49e28            # m^-3
e_C      = 1.602e-19          # C
c_light  = 2.998e8            # m/s
sigma_RT = 5.96e7             # S/m  (universal, 273 K)

def sigma_dc(T_K, RRR):
    """DC conductivity via Matthiessen's rule.

    Below 50 K phonons are frozen out and sigma_DC = RRR * sigma_RT — the
    simple 1/T phonon term is invalid there (umklapp/Bloch-Gruneisen) and
    must not be applied. Matches BuildDrudeMaterial in BBRMaterials.hh.
    """
    sig_imp = RRR * sigma_RT
    if T_K >= 50.0:
        sig_ph = sigma_RT * 273.0 / T_K
        return 1.0 / (1.0/sig_imp + 1.0/sig_ph)
    return sig_imp

def drude_tau(sig_dc_val):
    return sig_dc_val * m_e / (n_e * e_C**2)

def drude_D(freq_Hz, sig_dc_val):
    """Absorptance D = 1-R from the full complex Drude model.

    σ(ω) = σ_DC/(1−iωτ) kept complex: ε̃ = 1 + iσ/(ε₀ω), ñ = √ε̃,
    R = |(ñ−1)/(ñ+1)|². In the relaxation regime (ωτ >> 1, ω << ωp) this
    gives the flat absorptance D ≈ 2/(ωp·τ).
    """
    tau   = drude_tau(sig_dc_val)
    omega = 2.0 * np.pi * freq_Hz
    sigma_ac = sig_dc_val / (1.0 - 1j * omega * tau)
    eps_t = 1.0 + 1j * sigma_ac / (eps0 * omega)
    n_t   = np.sqrt(eps_t)
    R     = np.abs((n_t - 1.0) / (n_t + 1.0)) ** 2
    return 1.0 - np.clip(R, 0.0, 1.0)

def hagen_rubens_D(freq_Hz, sig_dc_val):
    omega = 2.0 * np.pi * freq_Hz
    return np.clip(2.0 * np.sqrt(2.0 * eps0 * omega / sig_dc_val), 0.0, 1.0)

out_dir = os.path.dirname(os.path.abspath(__file__))

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Drude scaling (two panels)
# ─────────────────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))
fig.subplots_adjust(wspace=0.35)

# ── Panel 1: D/D_HR vs ωτ ────────────────────────────────────────────────────
x = np.logspace(-2, 2, 1000)   # ωτ

# Exact result for the good-conductor limit (σ_DC/ε₀ω >> 1, our whole range)
# with the FULL complex σ(ω):
#   D/D_HR = √2 · cos(π/4 + arctan(ωτ)/2) · (1 + ω²τ²)^{1/4}
# → 1 at ωτ << 1 (Hagen-Rubens recovered)
# → 1/√(2ωτ) at ωτ >> 1 (relaxation regime, D → 2/(ωp·τ) flat in ω)
# H-R therefore OVERESTIMATES absorptance above f_break.
D_ratio = np.sqrt(2.0) * np.cos(np.pi/4.0 + np.arctan(x)/2.0) * (1.0 + x**2)**0.25
# Asymptotes
D_asym_lo = np.ones_like(x)
D_asym_hi = 1.0 / np.sqrt(2.0 * x)

ax1.loglog(x, D_ratio, 'k-', lw=2.2, zorder=4,
           label=r'Exact: $\sqrt{2}\cos\!\left(\frac{\pi}{4}+\frac{\arctan\omega\tau}{2}\right)(1+\omega^2\tau^2)^{1/4}$')
ax1.loglog(x[x < 1.5], D_asym_lo[x < 1.5], 'b--', lw=1.5, zorder=3,
           label=r'Low-$f$: $D = D_{\rm HR}$')
ax1.loglog(x[x > 0.7], D_asym_hi[x > 0.7], 'r--', lw=1.5, zorder=3,
           label=r'High-$f$: $1/\sqrt{2\omega\tau}$  ($D \to 2/\omega_p\tau$)')

ax1.axvline(1.0, color='gray', lw=1.0, ls='--', alpha=0.7)
ax1.annotate(r'$\omega\tau = 1$' + '\n' + r'$(f = f_{\rm break})$',
             xy=(1.0, 0.5), xytext=(2.5, 0.7), fontsize=8,
             arrowprops=dict(arrowstyle='->', color='gray', lw=0.8),
             color='gray', ha='left')

# Slope annotations
ax1.annotate(r'slope $0$', xy=(0.05, 1.05), fontsize=8,
             color='blue', ha='left')
ax1.annotate(r'slope $-1/2$', xy=(10, 0.3), fontsize=8,
             color='red', ha='left')

ax1.set_xlim(0.01, 100)
ax1.set_ylim(0.05, 2.0)
ax1.set_xlabel(r'$\omega\tau$')
ax1.set_ylabel(r'$D\,/\,D_{\rm HR}(\sigma_{\rm DC})$')
ax1.set_title('Universal Drude absorptance ratio')
ax1.legend(loc='lower left', fontsize=8, framealpha=0.9)
ax1.grid(True, which='both', ls=':', alpha=0.3)

# ── Panel 2: Re[σ] and Im[σ] vs ωτ ──────────────────────────────────────────
sigma_r_n = 1.0 / (1.0 + x**2)    # Re[σ(ω)]/σ_DC
sigma_i_n = x    / (1.0 + x**2)   # Im[σ(ω)]/σ_DC

ax2.loglog(x, sigma_r_n, color='steelblue', lw=2.2,
           label=r'${\rm Re}[\sigma(\omega)]/\sigma_{\rm DC} = 1/(1+\omega^2\tau^2)$')
ax2.loglog(x, sigma_i_n, color='darkorange', lw=2.2,
           label=r'${\rm Im}[\sigma(\omega)]/\sigma_{\rm DC} = \omega\tau/(1+\omega^2\tau^2)$')

# Asymptote guides for Re[σ]
ax2.loglog([0.01, 1.0], [1.0, 1.0],            'b:', lw=1.0, alpha=0.55)   # → 1
ax2.loglog([1.0, 100.0], [1.0, 1e-4],          'b:', lw=1.0, alpha=0.55)   # → (ωτ)^{-2}
# Asymptotes for Im[σ]
ax2.loglog([0.01, 1.0],  [0.01, 1.0],          color='darkorange', ls=':', lw=1.0, alpha=0.55)  # → ωτ
ax2.loglog([1.0, 100.0], [1.0, 1.0/100.0],     color='darkorange', ls=':', lw=1.0, alpha=0.55)  # → (ωτ)^{-1}

ax2.scatter([1.0], [0.5], s=50, color='darkorange', zorder=6,
            label=r'Peak of ${\rm Im}[\sigma]$ at $\omega\tau=1$: value $\sigma_{\rm DC}/2$')
ax2.axvline(1.0, color='gray', lw=1.0, ls='--', alpha=0.7)

# Slope annotations
ax2.annotate(r'slope $-2$', xy=(5.0, 0.02), fontsize=8, color='steelblue')
ax2.annotate(r'slope $+1$', xy=(0.02, 0.025), fontsize=8, color='darkorange')
ax2.annotate(r'slope $-1$', xy=(5.0, 0.12), fontsize=8, color='darkorange')

ax2.set_xlim(0.01, 100)
ax2.set_ylim(5e-5, 2.0)
ax2.set_xlabel(r'$\omega\tau$')
ax2.set_ylabel(r'$\sigma(\omega)\,/\,\sigma_{\rm DC}$')
ax2.set_title('Drude AC conductivity')
ax2.legend(loc='lower left', fontsize=7.5, framealpha=0.9)
ax2.grid(True, which='both', ls=':', alpha=0.3)

out1 = os.path.join(out_dir, "fig_drude_scaling.pdf")
fig.savefig(out1, bbox_inches='tight')
print(f"Saved: {out1}")
plt.close()

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Temperature dependence of D vs frequency
# ─────────────────────────────────────────────────────────────────────────────
fig2, ax = plt.subplots(1, 1, figsize=(6.5, 4.5))

freq_GHz = np.logspace(-1, 5, 700)
freq_Hz  = freq_GHz * 1e9

temps   = [4, 20, 77, 300]
cmap    = plt.cm.plasma
tcolors = [cmap(0.05), cmap(0.30), cmap(0.60), cmap(0.88)]

for T, col in zip(temps, tcolors):
    sig  = sigma_dc(T, RRR=100)
    tau  = drude_tau(sig)
    f_br = 1.0 / (2.0 * np.pi * tau) / 1e9
    D    = drude_D(freq_Hz, sig)
    lbl  = rf'$T = {T}\,\mathrm{{K}}$, $\sigma_{{DC}} = {sig:.1e}$ S/m, $f_{{\rm br}} = {f_br:.0f}$ GHz'
    ax.loglog(freq_GHz, D, color=col, lw=2.0, label=lbl)
    ax.axvline(f_br, color=col, lw=0.8, ls=':', alpha=0.6)

# H-R reference for 4 K
sig4  = sigma_dc(4, 100)
D_HR4 = hagen_rubens_D(freq_Hz, sig4)
ax.loglog(freq_GHz, D_HR4, 'k--', lw=1.2, alpha=0.55,
          label=r'H-R limit ($T=4\,\mathrm{K}$): $D \propto \omega^{1/2}$')

# Relaxation-regime plateau for 4 K: D → 2/(ωp·τ), flat in frequency
tau4   = drude_tau(sig4)
omega_p = np.sqrt(n_e * e_C**2 / (eps0 * m_e))
D_plateau = 2.0 / (omega_p * tau4)
plateau_str = f"{D_plateau:.1e}".replace("e-0", r"\times 10^{-") + "}"
ax.axhline(D_plateau, color='gray', lw=1.0, ls='-.',
           label=rf'Relaxation plateau ($T=4$ K): $D = 2/\omega_p\tau = {plateau_str}$')

ax.axvspan(50, 20000, alpha=0.07, color='gray', label='BBRsim range (50 GHz– 20 THz)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1, 1e5)
ax.set_ylim(2e-6, 0.12)
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel(r'Absorptance $D = 1 - R$')
ax.set_title('OFHC Cu (RRR = 100): Drude model, temperature dependence')
ax.set_xticks([1, 10, 100, 1000, 10000, 100000])
ax.set_xticklabels(['1 GHz', '10', '100', '1 THz', '10', '100 THz'])
ax.legend(fontsize=7.5, loc='upper left', framealpha=0.9)
ax.grid(True, which='both', ls=':', alpha=0.3)

plt.tight_layout()
out2 = os.path.join(out_dir, "fig_temp_dependence.pdf")
fig2.savefig(out2, bbox_inches='tight')
print(f"Saved: {out2}")
plt.close()
