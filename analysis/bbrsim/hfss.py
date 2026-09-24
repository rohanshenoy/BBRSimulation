"""bbrsim.hfss — Python mirror of the C++ HFSS diffraction sampler.

Reproduces, from the same CSV files, what ``BBRHFSSData`` and
``BBSimOpBoundaryProcess::HandleDiffractionBoundary`` compute for a photon
entering a ``vacuum_wg`` crack: the incidence-angle fold into the HFSS quarter
wedge, the polarization decomposition, the transmittance rule and the
outgoing-direction probability table. Validators use it to derive
deterministic predictions instead of running a second Monte Carlo.

Conventions (standard crack geometry; crack-local axes = world axes):
    normal_hat = +x   exit-face normal / propagation direction
    theta_hat  = +y   long dimension of the gap
    phi_hat    = +z   gap (short) dimension
HFSS incoming frame: z_i = normal_hat, x_i = phi_hat, y_i = -theta_hat.
    IWaveTheta = acos(-k . normal)                (180 deg = normal incidence)
    IWavePhi   = atan2(-k . theta_hat, k . phi_hat), folded into [0, 90] deg
                 by mirroring theta_hat (sign sy) and phi_hat (sign sx).
Outgoing frame (far_field.csv Theta, Phi), with theta_f = sy*theta_hat and
phi_f = sx*phi_hat:
    dir_out = sinT cosP normal + sinT sinP theta_f + cosT phi_f

Every formula here has a named counterpart in src/BBSimOpBoundaryProcess.cc
or src/BBRHFSSData.cc; keep them in lock-step.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

K_MIN_NORMAL_COMPONENT = 1e-6   # kMinNormalComponent in BBRHFSSData.cc
KEY_ROUND = 2                   # RoundDeg: incidence keys rounded to 0.01 deg


def default_base_dir():
    """Directory holding the ``<dataset>_Ephi=N`` folders: ``data/waveguides``
    of the checkout this file lives in.
    """
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(repo, "data", "waveguides")


@dataclass
class AngleDataset:
    """One (IWavePhi, IWaveTheta) incidence key of one crack dataset."""
    key: tuple                # (IWavePhi_deg, IWaveTheta_deg)
    T0: float                 # transmittance, Ephi=0 (theta-polarised) input, capped at 1
    T1: float                 # transmittance, Ephi=1 (phi-polarised) input, capped at 1
    theta_deg: np.ndarray     # far-field grid: outgoing Theta
    phi_deg: np.ndarray       # far-field grid: outgoing Phi
    F0_theta: np.ndarray      # complex rEtheta for Ephi=0 input
    F0_phi: np.ndarray        # complex rEphi   for Ephi=0 input
    F1_theta: np.ndarray      # complex rEtheta for Ephi=1 input
    F1_phi: np.ndarray        # complex rEphi   for Ephi=1 input


def load_dataset(dir_stem, base_dir=None):
    """Load ``<base>/<dir_stem>_Ephi={0,1}/{far_field,waveguide}.csv``.

    ``dir_stem`` is the directory name without the ``_Ephi=N`` suffix, i.e.
    ``<id>_<freq>GHz`` (e.g. ``InfParallelPlate_crack1Rohan_500GHz``) — the
    same string the C++ passes to the BBRHFSSData constructor.

    Returns {(IWavePhi, IWaveTheta): AngleDataset}. Mirrors BBRHFSSData's
    constructor: the Ephi=1 far-field rows are matched to the Ephi=0 rows by
    order within each incidence key, and raw power ratios above 1 are capped.
    """
    dataset_id = dir_stem   # local alias: the error messages below name the stem
    base = base_dir or default_base_dir()
    ff, wg = {}, {}
    for e in (0, 1):
        d = os.path.join(base, f"{dir_stem}_Ephi={e}")
        ff[e] = pd.read_csv(os.path.join(d, "far_field.csv"))
        wg[e] = pd.read_csv(os.path.join(d, "waveguide.csv"))

    out = {}
    keys = ff[0][["IWavePhi", "IWaveTheta"]].drop_duplicates()
    for phi_i, theta_i in keys.itertuples(index=False):
        key = (round(float(phi_i), KEY_ROUND), round(float(theta_i), KEY_ROUND))
        a = ff[0][(ff[0].IWavePhi == phi_i) & (ff[0].IWaveTheta == theta_i)].reset_index(drop=True)
        b = ff[1][(ff[1].IWavePhi == phi_i) & (ff[1].IWaveTheta == theta_i)].reset_index(drop=True)
        if len(a) != len(b) or not (np.allclose(a.Phi, b.Phi) and np.allclose(a.Theta, b.Theta)):
            raise ValueError(f"{dataset_id} key {key}: Ephi=0 and Ephi=1 far-field grids differ")
        T = {}
        for e in (0, 1):
            r = wg[e][(wg[e].IWavePhi == phi_i) & (wg[e].IWaveTheta == theta_i)].iloc[0]
            t = float(r.OutgoingPower) / float(r.IngoingPower) if r.IngoingPower > 0 else 0.0
            T[e] = min(1.0, t)   # load-time cap (BBRHFSSData::LoadWaveguide)
        out[key] = AngleDataset(
            key, T[0], T[1],
            a.Theta.to_numpy(float), a.Phi.to_numpy(float),
            a.rEtheta_real.to_numpy(float) + 1j * a.rEtheta_imag.to_numpy(float),
            a.rEphi_real.to_numpy(float) + 1j * a.rEphi_imag.to_numpy(float),
            b.rEtheta_real.to_numpy(float) + 1j * b.rEtheta_imag.to_numpy(float),
            b.rEphi_real.to_numpy(float) + 1j * b.rEphi_imag.to_numpy(float),
        )
    return out


def nearest_key(datasets, phi_deg, theta_deg):
    """Nearest incidence key by L2 distance in degrees (BBRHFSSData::FindDataset)."""
    return min(datasets, key=lambda k: (k[0] - phi_deg) ** 2 + (k[1] - theta_deg) ** 2)


@dataclass
class Incidence:
    """Folded incidence angles and the axes the sampler uses for this photon."""
    theta_deg: float
    phi_deg: float           # folded into [0, 90]
    phi_raw_deg: float       # before folding, in (-180, 180]
    sy: float                # -1 if theta_hat was mirrored
    sx: float                # -1 if phi_hat was mirrored
    normal: np.ndarray       # exit normal oriented along the photon
    theta_f: np.ndarray      # sy * theta_hat
    phi_f: np.ndarray        # sx * phi_hat


def fold_incidence(k, normal=(1., 0., 0.), theta_hat=(0., 1., 0.), phi_hat=(0., 0., 1.)):
    """Incoming angles + quarter-symmetry fold (HandleDiffractionBoundary).

    Transverse components below 1e-12 are snapped to +0 (as in the C++), so a
    photon with k_y == 0 exactly always folds with phi_raw = +180, not -180.
    """
    k = np.asarray(k, float)
    k = k / np.linalg.norm(k)
    n = np.asarray(normal, float)
    th = np.asarray(theta_hat, float)
    ph = np.asarray(phi_hat, float)
    if np.dot(k, n) < 0.:
        n = -n
    cos_val = min(1., max(-1., -np.dot(k, n)))
    theta_deg = float(np.degrees(np.arccos(cos_val)))
    # Snap tiny transverse components to +0 exactly as the C++ does, so a k with
    # k_y == 0 folds the same way regardless of signed-zero arithmetic.
    kt = float(-np.dot(k, th))
    kp = float(np.dot(k, ph))
    if abs(kt) < 1e-12:
        kt = 0.
    if abs(kp) < 1e-12:
        kp = 0.
    phi_raw = float(np.degrees(np.arctan2(kt, kp)))
    sy = -1. if phi_raw < 0. else 1.
    phi_deg = abs(phi_raw)
    sx = -1. if phi_deg > 90. else 1.
    if phi_deg > 90.:
        phi_deg = 180. - phi_deg
    return Incidence(theta_deg, phi_deg, phi_raw, sy, sx, n, sy * th, sx * ph)


def incoming_basis(inc):
    """(e_theta_in, e_phi_in): HFSS incoming spherical basis in world frame."""
    th, ph = np.radians(inc.theta_deg), np.radians(inc.phi_deg)
    e_theta = (np.sin(th) * inc.normal
               - np.cos(th) * np.sin(ph) * inc.theta_f
               + np.cos(th) * np.cos(ph) * inc.phi_f)
    e_phi = -np.cos(ph) * inc.theta_f - np.sin(ph) * inc.phi_f
    return e_theta, e_phi


def polarization_components(pol, inc):
    """(E_theta, E_phi) of a real polarization vector, normalised as in C++."""
    e_theta, e_phi = incoming_basis(inc)
    Et, Ep = float(np.dot(pol, e_theta)), float(np.dot(pol, e_phi))
    norm = float(np.hypot(Et, Ep))
    if norm > 1e-9:
        return Et / norm, Ep / norm
    return 1. / np.sqrt(2.), 1. / np.sqrt(2.)


def transmittance(ds, E_theta, E_phi):
    """Wang eq. 54 with the [0, 1] clamp (BBRHFSSData::GetTransmittance)."""
    return min(1., max(0., E_theta ** 2 * ds.T0 + E_phi ** 2 * ds.T1))


def direction_weights(ds, E_theta, E_phi):
    """Normalised sampling weights over the far-field grid rows.

    weight_i = sinT_i * |E_theta F0_i + E_phi F1_i|^2, zero for rows whose
    direction lies in the exit-face plane (BBRHFSSData::SampleOutgoingDirection).
    """
    sinT = np.sin(np.radians(ds.theta_deg))
    cosN = sinT * np.cos(np.radians(ds.phi_deg))
    Fth = E_theta * ds.F0_theta + E_phi * ds.F1_theta
    Fph = E_theta * ds.F0_phi + E_phi * ds.F1_phi
    w = sinT * (np.abs(Fth) ** 2 + np.abs(Fph) ** 2)
    w = np.where(cosN < K_MIN_NORMAL_COMPONENT, 0., w)
    s = w.sum()
    return w / s if s > 0 else np.full_like(w, 1. / len(w))


def outgoing_directions(ds, inc):
    """World-frame unit vectors of every far-field grid row, (N, 3)."""
    T = np.radians(ds.theta_deg)
    P = np.radians(ds.phi_deg)
    sinT, cosT, sinP, cosP = np.sin(T), np.cos(T), np.sin(P), np.cos(P)
    return (np.outer(sinT * cosP, inc.normal)
            + np.outer(sinT * sinP, inc.theta_f)
            + np.outer(cosT, inc.phi_f))


def random_polarization_mixture(ds, n_psi=360):
    """Predictions for uniformly random linear polarization.

    The gun draws (E_theta, E_phi) = (cos psi, sin psi) with psi uniform. A
    photon of polarization psi transmits with probability T(psi) and then
    samples its exit direction from its own normalised weight table, so the
    transmitted-photon direction distribution is the T-weighted average of the
    per-psi tables. Returns (mean transmittance, weights over grid rows).
    """
    psi = (np.arange(n_psi) + 0.5) / n_psi * 2. * np.pi
    W = np.zeros(len(ds.theta_deg))
    Tsum = 0.
    for p in psi:
        Et, Ep = float(np.cos(p)), float(np.sin(p))
        T = transmittance(ds, Et, Ep)
        if T <= 0.:
            continue
        W += T * direction_weights(ds, Et, Ep)
        Tsum += T
    T_mean = Tsum / n_psi
    return T_mean, (W / Tsum if Tsum > 0 else W)


def reflected_direction(k, normal=(1., 0., 0.)):
    """Specular reflection about the crack normal (the non-transmitted branch)."""
    k = np.asarray(k, float)
    n = np.asarray(normal, float)
    if np.dot(k, n) < 0.:
        n = -n
    return k - 2. * np.dot(k, n) * n


def binned_expectation(values, weights, edges):
    """Sum of ``weights`` falling in each bin of ``edges`` (for chi^2 tests)."""
    idx = np.clip(np.searchsorted(edges, values, side="right") - 1, 0, len(edges) - 2)
    out = np.zeros(len(edges) - 1)
    np.add.at(out, idx, weights)
    return out
