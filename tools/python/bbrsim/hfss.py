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
HFSS global frame in these axes: X = -phi_hat, Y = +theta_hat, Z = +normal_hat.
HFSS incidence angles are the spherical angles of the arrival direction r = -k:
    IWaveTheta = acos(-k . normal)                (180 deg = normal incidence)
    IWavePhi   = atan2(-k . theta_hat, k . phi_hat), folded into [0, 90] deg
                 by mirroring theta_hat (sign sy) and phi_hat (sign sx).
The exit coordinate system (x_e, y_e, z_e) = (Z, Y, -X) equals (normal,
theta_hat, phi_hat), so far_field.csv (Theta, Phi) and waveguide.csv (Y, Z) need
no transform. With theta_f = sy*theta_hat and phi_f = sx*phi_hat:
    dir_out = sinT cosP normal + sinT sinP theta_f + cosT phi_f

Every formula here has a named counterpart in library/src/BBSimOpBoundaryProcess.cc
or library/src/BBRHFSSData.cc; keep them in lock-step.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import paths

K_MIN_NORMAL_COMPONENT = 1e-6   # kMinNormalComponent in BBRHFSSData.cc
KEY_ROUND = 2                   # RoundDeg: incidence keys rounded to 0.01 deg


def default_base_dir():
    """<data root>/waveguides, the root resolved by bbrsim.paths.data_dir()."""
    return os.path.join(paths.data_dir(), "waveguides")


@dataclass
class AngleDataset:
    """One (IWavePhi, IWaveTheta) incidence key of one crack dataset."""
    key: tuple                # (IWavePhi_deg, IWaveTheta_deg)
    T0: float                 # transmittance, Ephi=0 (theta-polarised) input (normalized at load)
    T1: float                 # transmittance, Ephi=1 (phi-polarised) input (normalized at load)
    theta_deg: np.ndarray     # far-field grid: outgoing Theta
    phi_deg: np.ndarray       # far-field grid: outgoing Phi
    F0_theta: np.ndarray      # complex rEtheta for Ephi=0 input
    F0_phi: np.ndarray        # complex rEphi   for Ephi=0 input
    F1_theta: np.ndarray      # complex rEtheta for Ephi=1 input
    F1_phi: np.ndarray        # complex rEphi   for Ephi=1 input
    rho: complex = 0j         # normalized exit-field overlap <E0,E1> (polarization cross term)
    exit_y_m: np.ndarray = None   # waveguide.csv Y of every exit point (exit CS, metres; along theta_f)
    exit_z_m: np.ndarray = None   # waveguide.csv Z of every exit point (exit CS, metres; along phi_f)
    E0: np.ndarray = None         # (M, 3) complex exit field (Ex, Ey, Ez) for Ephi=0 input
    E1: np.ndarray = None         # (M, 3) complex exit field for Ephi=1 input, paired by row with E0


def load_dataset(dir_stem, base_dir=None):
    """Load ``<base>/<dir_stem>_Ephi={0,1}/{far_field,waveguide}.csv``.

    ``dir_stem`` is the directory name without the ``_Ephi=N`` suffix, i.e.
    ``<id>_<freq>GHz`` (e.g. ``InfParallelPlate_crack1Rohan_500GHz``) — the
    same string the C++ passes to the BBRHFSSData constructor.

    Returns {(IWavePhi, IWaveTheta): AngleDataset}. Mirrors BBRHFSSData's
    constructor: every row's Freq must agree with the stem's frequency to 0.1 %
    (BBR009); every numeric field must be a finite number (BBR013); the Ephi=1
    rows are matched to the Ephi=0 rows by order within each incidence key and
    must sit at the same far-field (Phi, Theta) and exit-point (X, Y, Z), with the
    same row count and no key of their own (BBR012); the two CSVs must have the
    same incidence keys (BBR007; BBR000 when neither has one); a key whose largest
    transmittance over linear polarizations exceeds 1 has T0 and T1 divided by it.
    The C++ errors are raised here as ValueError naming the code. No sidecar is
    read: the Python-side sidecar checks (F1-F13, C1-C5) are bbrsim.sidecar and
    validation/check_dataset_sidecars.py.
    """
    dataset_id = dir_stem   # local alias: the error messages below name the stem
    base = base_dir or default_base_dir()
    stem_token = dir_stem.rsplit("_", 1)[-1]
    f_stem = _strtod_full(stem_token[:-3]) if stem_token.endswith("GHz") else None
    if f_stem is None or not f_stem > 0.:
        raise ValueError(f"{dir_stem}: the stem does not end in _<freq>GHz")
    ff, wg = {}, {}
    for e in (0, 1):
        d = os.path.join(base, f"{dir_stem}_Ephi={e}")
        ff[e] = pd.read_csv(os.path.join(d, "far_field.csv"))
        wg[e] = pd.read_csv(os.path.join(d, "waveguide.csv"))
        for name, df in (("far_field.csv", ff[e]), ("waveguide.csv", wg[e])):
            for tok in df["Freq"].astype(str).unique():   # every row (BBR009)
                f = parse_frequency_GHz(tok)
                if f is None or abs(f - f_stem) > 1e-3 * f_stem:
                    raise ValueError(f"BBR009: Freq column {tok!r} in {d}/{name} does not match "
                                     f"the directory frequency {f_stem:g} GHz (tolerance 0.1 %)")
            # the fields the C++ parses with Num (all but Freq and Ephi) (BBR013)
            vals = df.drop(columns=["Freq", "Ephi"]).apply(pd.to_numeric, errors="coerce").to_numpy(float)
            if not np.isfinite(vals).all():
                raise ValueError(f"BBR013: non-numeric or non-finite field in {d}/{name}")

    # Incidence keys, in the order the C++ checks them: an Ephi=1 key without
    # Ephi=0 rows (BBR012, in the loaders), no key at all (BBR000), a key in one
    # CSV only (BBR007).
    kf, kw = _keys(ff[0]), _keys(wg[0])
    for name, df, k0 in (("far_field.csv", ff[1], kf), ("waveguide.csv", wg[1], kw)):
        extra = _keys(df) - k0
        if extra:
            raise ValueError(f"BBR012: {dataset_id}: Ephi=1 {name} has keys without Ephi=0 rows {sorted(extra)}")
    if not kf and not kw:
        raise ValueError(f"BBR000: {dataset_id}: no incidence key in either CSV")
    if kf != kw:
        raise ValueError(f"BBR007: {dataset_id}: far_field.csv and waveguide.csv have different "
                         f"incidence keys {sorted(kf ^ kw)}")

    out = {}
    keys = ff[0][["IWavePhi", "IWaveTheta"]].drop_duplicates()
    for phi_i, theta_i in keys.itertuples(index=False):
        key = (round(float(phi_i), KEY_ROUND), round(float(theta_i), KEY_ROUND))
        a = ff[0][(ff[0].IWavePhi == phi_i) & (ff[0].IWaveTheta == theta_i)].reset_index(drop=True)
        b = ff[1][(ff[1].IWavePhi == phi_i) & (ff[1].IWaveTheta == theta_i)].reset_index(drop=True)
        if len(a) != len(b) or not _same(a, b, ("Phi", "Theta"), 1e-9):
            raise ValueError(f"BBR012: {dataset_id} key {key}: Ephi=0 and Ephi=1 far-field grids differ")
        w = {e: wg[e][(wg[e].IWavePhi == phi_i) & (wg[e].IWaveTheta == theta_i)].reset_index(drop=True)
             for e in (0, 1)}
        if len(w[0]) != len(w[1]) or not _same(w[0], w[1], ("X", "Y", "Z"), 1e-12):
            raise ValueError(f"BBR012: {dataset_id} key {key}: Ephi=0 and Ephi=1 exit points differ")
        T, E = {}, {}
        for e in (0, 1):
            r = w[e].iloc[0]
            T[e] = float(r.OutgoingPower) / float(r.IngoingPower) if r.IngoingPower > 0 else 0.0
            E[e] = np.stack([w[e][f"E{c}_real"].to_numpy(float) + 1j * w[e][f"E{c}_imag"].to_numpy(float)
                             for c in "xyz"], axis=1)
        # Polarization cross term (BBRHFSSData constructor): the Ephi=1 exit points pair
        # with the Ephi=0 ones by order; rho = sum E0.E1* / sqrt(sum|E0|^2 sum|E1|^2).
        p0, p1 = float(np.sum(np.abs(E[0]) ** 2)), float(np.sum(np.abs(E[1]) ** 2))
        rho = complex(np.sum(E[0] * np.conj(E[1])) / np.sqrt(p0 * p1)) if p0 > 0 and p1 > 0 else 0j
        # Load-time normalization (BBRHFSSData constructor): the largest T over linear
        # polarizations, the top eigenvalue of [[T0, c], [c, T1]] with c = sqrt(T0 T1) Re rho,
        # above 1 (the HFSS port-normalization artefact) divides both T0 and T1.
        c = np.sqrt(T[0] * T[1]) * rho.real
        lam = float(0.5 * (T[0] + T[1]) + np.hypot(0.5 * (T[0] - T[1]), c))
        if lam > 1.0:
            T = {e: T[e] / lam for e in (0, 1)}
        out[key] = AngleDataset(
            key, T[0], T[1],
            a.Theta.to_numpy(float), a.Phi.to_numpy(float),
            a.rEtheta_real.to_numpy(float) + 1j * a.rEtheta_imag.to_numpy(float),
            a.rEphi_real.to_numpy(float) + 1j * a.rEphi_imag.to_numpy(float),
            b.rEtheta_real.to_numpy(float) + 1j * b.rEtheta_imag.to_numpy(float),
            b.rEphi_real.to_numpy(float) + 1j * b.rEphi_imag.to_numpy(float),
            rho,
            w[0].Y.to_numpy(float), w[0].Z.to_numpy(float), E[0], E[1],
        )
    return out


def _keys(df):
    """The incidence keys of one CSV, rounded as the C++ MakeKey rounds them."""
    return {(round(float(p), KEY_ROUND), round(float(t), KEY_ROUND))
            for p, t in df[["IWavePhi", "IWaveTheta"]].drop_duplicates().itertuples(index=False)}


def _same(a, b, cols, tol):
    """Row-by-row equality of two frames' columns to tol (the C++ BBR012 tolerances:
    1e-9 deg for Phi/Theta, 1e-12 m for X/Y/Z)."""
    return all(np.all(np.abs(a[c].to_numpy(float) - b[c].to_numpy(float)) <= tol) for c in cols)


def nearest_key(datasets, phi_deg, theta_deg):
    """Nearest incidence key by L2 distance in degrees (BBRHFSSData::FindDataset).

    The keys are scanned in ascending (IWavePhi, IWaveTheta) order, the std::map
    order of the C++, and the first minimum wins, so a tie resolves the same way.
    """
    return min(sorted(datasets), key=lambda k: (k[0] - phi_deg) ** 2 + (k[1] - theta_deg) ** 2)


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
    """(e_theta_in, e_phi_in) of HandleDiffractionBoundary in the world frame: exactly
    (-e_theta, -e_phi) of the HFSS spherical basis at the arrival direction r = -k, written in
    the folded frame (X, Y, Z) = (-phi_f, theta_f, normal). The common sign cancels in T, in the
    cross term and in |E|^2 and flips pol_out only, which is the same state; the +sin(theta)
    term is the -e_theta component, not a sign error."""
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
    """Wang eq. 58 applied to eq. 53 (BBRHFSSData::GetTransmittance): the transmitted power
    of a mixed polarization includes the cross term 2 Et Ep sqrt(T0 T1) Re rho."""
    T = (E_theta ** 2 * ds.T0 + E_phi ** 2 * ds.T1
         + 2. * E_theta * E_phi * np.sqrt(ds.T0 * ds.T1) * ds.rho.real)
    return min(1., max(0., T))


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


def exit_position_weights(ds, E_theta, E_phi):
    """Normalised sampling weights over the exit points (rows of waveguide.csv).

    weight_j = |E_theta E0_j + E_phi E1_j|^2 summed over (Ex, Ey, Ez), the runtime CDF
    of BBRHFSSData::SampleExitPosition; uniform when every weight is zero, as the C++
    then draws a uniform index. Point j sits at ds.exit_y_m[j] along theta_f and
    ds.exit_z_m[j] along phi_f from the exit-face centre.
    """
    w = np.sum(np.abs(E_theta * ds.E0 + E_phi * ds.E1) ** 2, axis=1)
    s = w.sum()
    return w / s if s > 0 else np.full(len(w), 1. / len(w))


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


# ---------------------------------------------------------------------------
# Frequency grid: discovery and selection
#
# Mirrors BBRCrackLibrary::Discover and the selection block of
# BBRCrackLibrary::Lookup. Keep the two implementations in lock-step: the
# validators compare the C++ choice, recorded per photon in the crossings
# column hfss_freq_GHz, against select_frequency() evaluated on the same
# directory tree.
# ---------------------------------------------------------------------------

H_EV_S = 4.135667696e-15   # Planck constant in eV s (CODATA 2018 = CLHEP's value)


def photon_frequency_GHz(energy_eV):
    """Photon frequency in GHz from its energy in eV (scalar or array)."""
    nu = np.asarray(energy_eV, dtype=float) / H_EV_S / 1e9
    return float(nu) if nu.ndim == 0 else nu


# What C strtod parses in full (C locale), as the C++ full-parse checks accept it
# (end == token end): optional leading whitespace, a sign, then a decimal number,
# a hexadecimal one (0x, optionally with a binary exponent p), inf, infinity or nan
# in any case, and nothing after it. Python's float() differs: it also takes
# trailing whitespace and digit underscores, and no hex.
_STRTOD_FULL = re.compile(
    r"[ \t\n\v\f\r]*[+-]?(?:"
    r"(?P<hex>0[xX](?:[0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)(?:[pP][+-]?[0-9]+)?)"
    r"|(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    r"|(?i:inf|infinity|nan))")


def _strtod_full(token):
    """float(token) when C strtod would parse all of it, else None."""
    m = _STRTOD_FULL.fullmatch(token)
    if m is None:
        return None
    t = token.lstrip(" \t\n\v\f\r")
    return float.fromhex(t) if m.group("hex") else float(t)


def parse_frequency_GHz(token):
    """BBRHFSSData::ParseFrequencyGHz: '<number>MHz|GHz|THz' -> GHz, or None.

    Whitespace anywhere is ignored and the unit is case-insensitive; the number
    must parse in full and be > 0.
    """
    t = re.sub(r"[ \t\n\v\f\r]", "", str(token))
    if len(t) < 4:
        return None
    mult = {"mhz": 1e-3, "ghz": 1., "thz": 1e3}.get(t[-3:].lower())
    v = _strtod_full(t[:-3]) if mult is not None else None
    return v * mult if v is not None and v > 0. else None


def _parse_stem_token(name, prefix, suffix):
    """(value, token) of '<prefix><token><suffix>', or None if it does not match.

    The token must parse in full as C strtod does and be > 0, as in
    BBRCrackLibrary::Discover (so ' 500', '0x1f4' and 'inf' count; '500 ',
    '5_00', 'nan', '0' and '-5' do not).
    """
    if not (name.startswith(prefix) and name.endswith(suffix)):
        return None
    token = name[len(prefix):len(name) - len(suffix)]
    value = _strtod_full(token)
    return (value, token) if value is not None and value > 0. else None


def discover_frequencies(dataset_id, base_dir=None):
    """[(freq_GHz, dir_stem), ...] for '<dataset_id>_<freq>GHz_Ephi=0' folders.

    Sorted by frequency. The token is kept verbatim because a parsed double
    cannot regenerate it ('1.5e3' and '1500' are the same value but different
    directory names). Raises on a legacy '<id>_Ephi=0' folder, on two folders
    whose tokens parse to the same value, and when nothing matches.
    """
    base = base_dir or default_base_dir()
    prefix, suffix = f"{dataset_id}_", "GHz_Ephi=0"
    found = {}
    for name in sorted(os.listdir(base)):
        if not os.path.isdir(os.path.join(base, name)):
            continue
        if name == f"{dataset_id}_Ephi=0":
            raise ValueError(
                f"legacy directory {name} has no frequency: "
                f"rename it to {prefix}<freq>{suffix}")
        parsed = _parse_stem_token(name, prefix, suffix)
        if parsed is None:
            continue
        value, token = parsed
        if value in found:
            raise ValueError(
                f"duplicate frequency {value} GHz for {dataset_id}: "
                f"{found[value]} and {dataset_id}_{token}GHz")
        found[value] = f"{dataset_id}_{token}GHz"
    if not found:
        raise FileNotFoundError(
            f"no {prefix}<freq>{suffix} directories under {base}")
    return sorted(found.items())


def select_frequency(entries, nu_GHz):
    """(freq_GHz, dir_stem, clamped) for a photon of frequency ``nu_GHz``.

    ``entries`` is the list returned by discover_frequencies (ascending).
    clamped is -1 below the grid, +1 above it, 0 on-grid. A single-entry grid
    is never clamped (there is nothing to choose), which keeps single-frequency
    real data silent. Ties in log distance go to the lower frequency.
    """
    if len(entries) == 1:
        f, stem = entries[0]
        return f, stem, 0
    if nu_GHz < entries[0][0]:
        f, stem = entries[0]
        return f, stem, -1
    if nu_GHz > entries[-1][0]:
        f, stem = entries[-1]
        return f, stem, +1
    lnu = np.log10(nu_GHz)
    best, k = float("inf"), 0
    for i, (f, _) in enumerate(entries):   # ascending; strict '<' keeps the lower on a tie
        d = abs(np.log10(f) - lnu)
        if d < best:
            best, k = d, i
    f, stem = entries[k]
    return f, stem, 0


def load_frequency_set(dataset_id, base_dir=None):
    """{freq_GHz: {(IWavePhi, IWaveTheta): AngleDataset}} for every frequency."""
    base = base_dir or default_base_dir()
    return {f: load_dataset(stem, base)
            for f, stem in discover_frequencies(dataset_id, base)}
