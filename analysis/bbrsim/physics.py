"""Shared BBRsim physics: Drude copper reflectance, Planck spectrum, Hagen-Rubens.

Single source of truth for the theory previously copy-pasted across the
scripts/ validators and plots. Mirrors the C++ BBRMaterials (Drude) and
ThermalSurface (Planck) math. Pure functions, NumPy-vectorized, no I/O.
"""
import numpy as np
from scipy.integrate import quad

# ── physical constants (SI unless noted) ─────────────────────────────────────
SIGMA_RT = 5.96e7             # S/m, universal Cu conductivity at 273 K
M_E      = 9.109e-31          # kg, electron mass
N_E      = 8.49e28            # m^-3, Cu free-electron density
E_CHARGE = 1.602e-19          # C
EPS0     = 8.8541878128e-12   # F/m
H_EVS    = 4.13566769692e-15  # eV*s, Planck constant
K_EV     = 8.617333e-5        # eV/K, Boltzmann constant
PLANCK_PEAK_U = 1.5936        # peak of u^2/(e^u-1), u = E/kT


def ev_to_hz(E_eV):
    """Photon energy [eV] -> frequency [Hz]."""
    return np.asarray(E_eV, dtype=float) / H_EVS


def hz_to_ev(freq_Hz):
    """Frequency [Hz] -> photon energy [eV]."""
    return np.asarray(freq_Hz, dtype=float) * H_EVS


def sigma_dc(RRR, T_K=4.0):
    """DC conductivity [S/m] via Matthiessen's rule.

    Below 50 K phonons are frozen: sigma = RRR * SIGMA_RT (impurity term).
    At/above 50 K add the linear-in-1/T phonon term sigma_ph = SIGMA_RT*273/T.
    """
    sigma_imp = RRR * SIGMA_RT
    if T_K >= 50.0:
        sigma_ph = SIGMA_RT * 273.0 / T_K
        return 1.0 / (1.0 / sigma_imp + 1.0 / sigma_ph)
    return sigma_imp


def drude_tau(RRR, T_K=4.0):
    """Drude relaxation time [s] from the DC conductivity."""
    return sigma_dc(RRR, T_K) * M_E / (N_E * E_CHARGE**2)


def drude_reflectance(freq_Hz, RRR, T_K=4.0):
    """Normal-incidence reflectance R from the full COMPLEX Drude model.

    sigma(w) = sigma_DC/(1 - i w tau) kept complex; eps = 1 + i sigma/(eps0 w);
    n = sqrt(eps); R = |(n-1)/(n+1)|^2. The real-sigma closed form is NOT used:
    it drops the plasma term in Re(eps) and overestimates D ~30x at 500 GHz
    (the bug fixed 2026-06-09). Vectorized over freq_Hz.
    """
    freq_Hz = np.asarray(freq_Hz, dtype=float)
    sdc   = sigma_dc(RRR, T_K)
    tau   = drude_tau(RRR, T_K)
    omega = 2.0 * np.pi * freq_Hz
    sigma = sdc / (1.0 - 1j * omega * tau)
    eps   = 1.0 + 1j * sigma / (EPS0 * omega)
    n     = np.sqrt(eps)
    R     = np.abs((n - 1.0) / (n + 1.0)) ** 2
    return np.clip(R, 0.0, 1.0)


def drude_absorptance(freq_Hz, RRR, T_K=4.0):
    """Absorptance D = 1 - R from the full complex Drude model."""
    return 1.0 - drude_reflectance(freq_Hz, RRR, T_K)


def hagen_rubens_absorptance(freq_Hz, sigma_SI):
    """Low-frequency-limit absorptance D = 2*sqrt(2 eps0 omega / sigma).

    Valid only for omega*tau << 1; overestimates D in the relaxation regime.
    """
    omega = 2.0 * np.pi * np.asarray(freq_Hz, dtype=float)
    return 2.0 * np.sqrt(2.0 * EPS0 * omega / sigma_SI)


def planck_photon_number_pdf(E_eV, T_K):
    """Unnormalized Planck photon-NUMBER spectral density E^2/(e^{E/kT}-1).

    This is the spectrum the emitter samples (photon number, not energy).
    Vectorized; returns 0 where the exponent would overflow.
    """
    E_eV = np.asarray(E_eV, dtype=float)
    x = E_eV / (K_EV * T_K)
    valid = (E_eV > 0) & (x < 500.0)
    safe = np.clip(x, 1e-300, 500.0)   # avoid expm1(0)=0 division on masked entries
    return np.where(valid, E_eV**2 / np.expm1(safe), 0.0)


def planck_weighted_absorptance(RRR, T_K=4.0, f_lo_Hz=10e9, f_hi_Hz=20e12):
    """Planck-photon-number-weighted Cu absorptance over [f_lo, f_hi].

    A = int D(f) pdf(E(f)) df / int pdf(E(f)) df, with E = h f. This is the
    average absorptance a Planck-spectrum photon beam sees on the Cu wall.
    """
    def num(f):
        return float(drude_absorptance(f, RRR, T_K)
                     * planck_photon_number_pdf(H_EVS * f, T_K))

    def den(f):
        return float(planck_photon_number_pdf(H_EVS * f, T_K))

    n, _ = quad(num, f_lo_Hz, f_hi_Hz, limit=200)
    d, _ = quad(den, f_lo_Hz, f_hi_Hz, limit=200)
    return n / d
