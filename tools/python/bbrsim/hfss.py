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

import math
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
    schema: int = 1               # schema major of the dataset's sidecar (1 without one): ExitPoint-less twin of fSchema
    w: np.ndarray = None          # (M,) exit quadrature weights, ExitPoint::w (all ones under 1.x)
    dark0: bool = False           # BBRKeyRecord::dark0: the Ephi=0 basis is a numerical zero (2.x, OR-038)
    dark1: bool = False           # BBRKeyRecord::dark1


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
    transmittance over linear polarizations exceeds 1 has T0 and T1 divided by it;
    a basis polarization with T > 0 needs outward far-field intensity and exit field
    (BBR027).
    The C++ errors are raised here as ValueError naming the code. The Python-side sidecar
    checks (F1-F13, C1-C5) are bbrsim.sidecar and validation/check_dataset_sidecars.py; this
    function reads ``<base>/<dir_stem>.dataset.json`` when it exists only to learn the schema
    major (``BBRHFSSData`` takes ``schemaMajor`` from the sidecar the same way) and, under 2.x,
    the exit lattice and the per-key dark flags.

    Schema 2.x (OR-034, OR-038, HF-040, HF-042; BBRHFSSData's ``fSchema == 2`` branches):
    waveguide.csv is read by the 2.0 positions, Ex = (9, 10), Ey = (11, 12), Ez = (13, 14) as
    (re, im) and InterceptFactor = 15; T_p = OutgoingPower / (IngoingPower * InterceptFactor) from the
    key's first row, kept above 1 with no eigenvalue rescale (``GetTransmittance`` clamps per photon);
    each exit point gets the trapezoid weight ``w`` of its lattice node (AssignQuadrature); rho is
    over (Ex, Ey) with those weights, 0 if either basis is dark. A 2.x CSV never falls back to
    the 1.0 column map: a 16-column header without a 2.x sidecar is BBR024.
    """
    dataset_id = dir_stem   # local alias: the error messages below name the stem
    base = base_dir or default_base_dir()
    stem_token = dir_stem.rsplit("_", 1)[-1]
    f_stem = _strtod_full(stem_token[:-3]) if stem_token.endswith("GHz") else None
    if f_stem is None or not f_stem > 0.:
        raise ValueError(f"{dir_stem}: the stem does not end in _<freq>GHz")
    from . import sidecar   # lazy: bbrsim.sidecar imports this module
    sc = sidecar.load(base, dir_stem) if os.path.exists(sidecar.path_for(base, dir_stem)) else None
    major = sidecar.schema_major(sc) if sc is not None else 1   # BBR024 for a version BBRsim does not read
    v2 = major == 2
    if v2:   # E_t below is the HFSS-global (Ex, Ey); the C++ refuses another frame at parse (CDX-035)
        sidecar.check_field_frame(sc, major)
    ff, wg = {}, {}
    for e in (0, 1):
        d = os.path.join(base, f"{dir_stem}_Ephi={e}")
        ff[e] = pd.read_csv(os.path.join(d, "far_field.csv"))
        wg[e] = pd.read_csv(os.path.join(d, "waveguide.csv"))
        wide = wg[e].shape[1] == len(sidecar.WG_COLUMNS_2)
        if wide and sc is None:   # never the 1.0 map on a 2.0 table: it would read Re Ex, Im Ex, Im Ey, Re Ez as E_t
            raise ValueError(f"BBR024: {d}/waveguide.csv has the 16-column schema-2 header but there is no "
                             f"{dir_stem}.dataset.json beside it; BBRsim reads a 2.x table only with its sidecar")
        if wide != v2 or (v2 and list(wg[e].columns) != sidecar.WG_COLUMNS_2):   # CheckHeader (C1)
            raise ValueError(f"BBR013: {d}/waveguide.csv header {list(wg[e].columns)} is not the "
                             f"{'2.0' if v2 else '1.0'} column list the sidecar declares")
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
        rec = _key_record(sc, (float(phi_i), float(theta_i))) if v2 else None
        T, E = {}, {}
        for e in (0, 1):
            r = w[e].iloc[0]
            if v2:
                # LoadWaveguide, 2.x: T = Out / (In * InterceptFactor) from the key's first row, kept above 1
                f = float(r.InterceptFactor)
                T[e] = float(r.OutgoingPower) / (float(r.IngoingPower) * f) if float(r.IngoingPower) * f > 0 else 0.0
                # the 2.0 positions Ex = (9, 10), Ey = (11, 12), Ez = (13, 14) as (re, im), by position
                col = w[e].iloc[:, 9:15].to_numpy(float)
                E[e] = np.stack([col[:, 0] + 1j * col[:, 1], col[:, 2] + 1j * col[:, 3],
                                 col[:, 4] + 1j * col[:, 5]], axis=1)
            else:
                T[e] = float(r.OutgoingPower) / float(r.IngoingPower) if r.IngoingPower > 0 else 0.0
                E[e] = np.stack([w[e][f"E{c}_real"].to_numpy(float) + 1j * w[e][f"E{c}_imag"].to_numpy(float)
                                 for c in "xyz"], axis=1)
        if v2:
            wq = _quadrature(sc, w[0], dataset_id, key)
            dark = (rec["dark0"], rec["dark1"])
            # Polarization cross term, 2.x (HF-042): rho over the tangential (Ex, Ey) with the quadrature
            # weights, the sampler's own density; 0 if either basis is dark. No rescale: values above 1 are
            # kept and transmittance() caps per photon (HF-039(2)).
            rho = 0j
            if not any(dark):
                ta, tb = E[0][:, :2], E[1][:, :2]
                p0, p1 = float(np.sum(wq[:, None] * np.abs(ta) ** 2)), float(np.sum(wq[:, None] * np.abs(tb) ** 2))
                if p0 > 0 and p1 > 0:
                    rho = complex(np.sum(wq[:, None] * ta * np.conj(tb)) / np.sqrt(p0 * p1))
        else:
            wq, dark = np.ones(len(w[0])), (False, False)
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
        # BBR027 (BBRHFSSData constructor): a transmitting basis polarization needs outward far-field
        # intensity and exit field somewhere (CDX-008)
        sinT = np.sin(np.radians(a.Theta.to_numpy(float)))
        outward = sinT * np.cos(np.radians(a.Phi.to_numpy(float))) >= K_MIN_NORMAL_COMPONENT
        for e, f in ((0, a), (1, b)):
            if T[e] > 0:
                ff_w = float(np.sum((sinT * (f.rEtheta_real ** 2 + f.rEtheta_imag ** 2 + f.rEphi_real ** 2
                                             + f.rEphi_imag ** 2).to_numpy(float))[outward]))
                # exit density: 1.x |E|^2 over (Ex, Ey, Ez); 2.x w |E_t|^2 over (Ex, Ey)
                ex_w = float(np.sum(wq[:, None] * np.abs(E[e][:, :2]) ** 2)) if v2 else float(np.sum(np.abs(E[e]) ** 2))
                if not ff_w > 0 or not ex_w > 0:
                    what = "exit field is zero at every exit point" if ff_w > 0 else "far field has no outward intensity"
                    raise ValueError(f"BBR027: {dataset_id} key {key}: T{e} = {T[e]!r} > 0, but its {what}")
        out[key] = AngleDataset(
            key, T[0], T[1],
            a.Theta.to_numpy(float), a.Phi.to_numpy(float),
            a.rEtheta_real.to_numpy(float) + 1j * a.rEtheta_imag.to_numpy(float),
            a.rEphi_real.to_numpy(float) + 1j * a.rEphi_imag.to_numpy(float),
            b.rEtheta_real.to_numpy(float) + 1j * b.rEtheta_imag.to_numpy(float),
            b.rEphi_real.to_numpy(float) + 1j * b.rEphi_imag.to_numpy(float),
            rho,
            w[0].Y.to_numpy(float), w[0].Z.to_numpy(float), E[0], E[1],
            major, wq, bool(dark[0]), bool(dark[1]),
        )
    return out


def _key_record(sc, key):
    """The 2.x per-key record of the CSV key ``key`` = (raw IWavePhi, raw IWaveTheta): matched on the raw
    angles at 1e-9 degree, BBRDatasetSidecar::Key's rule and the one LoadWaveguide applies per row, not on
    the 0.01 degree rounded group key (which C2 uses). A key without a record is BBR012."""
    from . import sidecar
    for r in sidecar.key_records(sc):
        if abs(r["phi"] - key[0]) <= sidecar.KEY_TOL and abs(r["theta"] - key[1]) <= sidecar.KEY_TOL:
            return r
    raise ValueError(f"BBR012: key ({key[0]!r}, {key[1]!r}) has no per-key record in the sidecar (angles matched "
                     "to 1e-9 degree against the declared excitation grid)")


def _axis_index(v, axis, name, key):
    """AxisNode (BBRHFSSData.cc): the index i of the declared node min + i (max - min)/(count - 1) that v
    matches to 1e-6 of the step, else BBR012."""
    n, lo, hi = int(axis["count"]), float(axis["min"]), float(axis["max"])
    if n < 2:
        raise ValueError(f"BBR024: key {key}: exit_field.grid has count {n} on the {name} axis; a rectangle axis "
                         "needs at least 2 nodes (count < 2 is refused at parse)")
    step = (hi - lo) / (n - 1)
    i = int(math.floor((v - lo) / step + 0.5))   # std::lround
    if 0 <= i < n and abs(v - (lo + i * step)) <= 1e-6 * step:
        return i, n
    raise ValueError(f"BBR012: key {key}: exit coordinate {name} = {v!r} m is not a node of the declared lattice; "
                     "quadrature weights are assigned on the unclamped lattice")


def _quadrature(sc, wg0, dataset_id, key):
    """ExitPoint::w of every exit point (BBRHFSSData::AssignQuadrature, HF-040): a rectangle's trapezoid by
    grid index, a factor 1/2 per axis on its first and last node (a corner 1/4); a disc ("uniform_disc") 1.
    Dimensionless, from the unclamped positions."""
    ef = sc["exit_field"]
    if ef["cross_section"]["shape"] != "rectangle":
        return np.ones(len(wg0))
    out = np.empty(len(wg0))
    for j, (y, z) in enumerate(zip(wg0.Y.to_numpy(float), wg0.Z.to_numpy(float))):
        wy = wz = 1.0
        for v, name, which in ((y, "Y", "y_e"), (z, "Z", "z_e")):
            i, n = _axis_index(v, ef["grid"][which], name, key)
            f = 0.5 if (i == 0 or i == n - 1) else 1.0
            if which == "y_e":
                wy = f
            else:
                wz = f
        out[j] = wy * wz
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
    of a mixed polarization includes the cross term 2 Et Ep sqrt(T0 T1) Re rho. Under 2.x T0 and T1
    are kept above 1 and this clamp is the per-photon cap (HF-039(2)); a dark basis has T = 0 and
    rho = 0, so it adds nothing (OR-038)."""
    T = (E_theta ** 2 * ds.T0 + E_phi ** 2 * ds.T1
         + 2. * E_theta * E_phi * np.sqrt(ds.T0 * ds.T1) * ds.rho.real)
    return min(1., max(0., T))


def direction_weights(ds, E_theta, E_phi):
    """Normalised sampling weights over the far-field grid rows.

    weight_i = sinT_i * |E_theta F0_i + E_phi F1_i|^2, zero for rows whose
    direction lies in the exit-face plane (BBRHFSSData::SampleOutgoingDirection). A dark basis
    (2.x, OR-038) enters no mixture: its stored far field is numerical residue.
    """
    cT = 0. if ds.dark0 else E_theta   # cT, cP of SampleOutgoingDirection
    cP = 0. if ds.dark1 else E_phi
    sinT = np.sin(np.radians(ds.theta_deg))
    cosN = sinT * np.cos(np.radians(ds.phi_deg))
    Fth = cT * ds.F0_theta + cP * ds.F1_theta
    Fph = cT * ds.F0_phi + cP * ds.F1_phi
    w = sinT * (np.abs(Fth) ** 2 + np.abs(Fph) ** 2)
    w = np.where(cosN < K_MIN_NORMAL_COMPONENT, 0., w)
    s = w.sum()
    if not s > 0:   # the basis far fields cancel: BBRHFSSData raises BBR027, it never draws uniformly
        raise ValueError(f"BBR027: no outward far-field intensity for (E_theta, E_phi) = ({E_theta}, {E_phi}) "
                         f"at key {ds.key}")
    return w / s


def exit_position_weights(ds, E_theta, E_phi):
    """Normalised sampling weights over the exit points (rows of waveguide.csv).

    weight_j = |E_theta E0_j + E_phi E1_j|^2 summed over (Ex, Ey, Ez), the runtime CDF
    of BBRHFSSData::SampleExitPosition; BBR027 when every weight is zero, as the C++
    then refuses to draw. Point j sits at ds.exit_y_m[j] along theta_f and
    ds.exit_z_m[j] along phi_f from the exit-face centre.

    Schema 2.x (HF-040, SampleExitPosition with fSchema == 2): weight_j = w_j |cT E_t0,j + cP E_t1,j|^2
    over the tangential (Ex, Ey) with the quadrature weight w_j; a dark basis enters no mixture
    (cT or cP = 0, OR-038).
    """
    cT = 0. if ds.dark0 else E_theta   # cT, cP of SampleExitPosition
    cP = 0. if ds.dark1 else E_phi
    if ds.schema == 2:
        w = ds.w * np.sum(np.abs(cT * ds.E0[:, :2] + cP * ds.E1[:, :2]) ** 2, axis=1)
    else:
        w = np.sum(np.abs(cT * ds.E0 + cP * ds.E1) ** 2, axis=1)
    s = w.sum()
    if not s > 0:
        raise ValueError(f"BBR027: no exit field for (E_theta, E_phi) = ({E_theta}, {E_phi}) at key {ds.key}")
    return w / s


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
    per-psi tables. Returns (mean transmittance, weights over grid rows). Each psi's T is the
    clamped one of GetTransmittance, so under 2.x, where T0 and T1 may exceed 1, the mean is the
    per-photon-capped one.
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

# CLHEP's units and Planck constant, built in CLHEP's own order
# (CLHEP/Units/SystemOfUnits.h, PhysicalConstants.h): electronvolt = 1e-6 MeV,
# joule = electronvolt / e_SI, second = 1e9 ns, h_Planck = 6.62607015e-34 * joule
# * second in MeV ns. photon_frequency_GHz divides by it as BBSimOpBoundaryProcess
# does, so a frequency computed here equals the C++ one bit for bit and the
# validators' selection checks are exact.
_E_SI = 1.602176634e-19                                # CLHEP::e_SI, coulomb
_ELECTRONVOLT = 1.e-6 * 1.                             # CLHEP::electronvolt, MeV
_JOULE = _ELECTRONVOLT / _E_SI                         # CLHEP::joule, MeV
_SECOND = 1.e+9 * 1.                                   # CLHEP::second, ns
H_PLANCK_MEV_NS = 6.62607015e-34 * _JOULE * _SECOND    # CLHEP::h_Planck, MeV ns
H_EV_S = H_PLANCK_MEV_NS / (_ELECTRONVOLT * _SECOND)   # the same in eV s, 4.135667696923858e-15


def photon_frequency_GHz(energy_eV):
    """Photon frequency in GHz from its energy in eV (scalar or array), computed as
    BBSimOpBoundaryProcess does: E[MeV] / h_Planck[MeV ns], where 1/ns is 1 GHz."""
    nu = (np.asarray(energy_eV, dtype=float) * _ELECTRONVOLT) / H_PLANCK_MEV_NS
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


_DATASET_DIR = re.compile(r"(?P<id>.+)_(?P<token>[^_]+)GHz_Ephi=0")


def discover_datasets(base_dir=None):
    """{dataset_id: [(freq_GHz, dir_stem), ...]} for every dataset under ``base_dir``.

    Each '<id>_<token>GHz_Ephi=0' folder whose token parses as BBRCrackLibrary
    parses it names a dataset; each ID's grid is then exactly what
    discover_frequencies returns, with its errors (two tokens of one value).
    A '<id>_Ephi=0' folder with no frequency raises, as for one ID. Folders with
    a token that does not parse, and files, are skipped. IDs come sorted.
    """
    base = base_dir or default_base_dir()
    ids = set()
    for name in sorted(os.listdir(base)):
        if not name.endswith("_Ephi=0") or not os.path.isdir(os.path.join(base, name)):
            continue
        if not name[:-len("_Ephi=0")].endswith("GHz"):
            raise ValueError(f"legacy directory {name} has no frequency: rename it to "
                             f"<id>_<freq>GHz_Ephi=0")
        m = _DATASET_DIR.fullmatch(name)
        if m and _parse_stem_token(name, f"{m.group('id')}_", "GHz_Ephi=0") is not None:
            ids.add(m.group("id"))
    return {i: discover_frequencies(i, base) for i in sorted(ids)}


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
