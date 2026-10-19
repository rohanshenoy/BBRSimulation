"""bbrsim.sidecar - the HFSS dataset sidecar <stem>.dataset.json (schema 1.x and 2.x).

Python twin of library/src/BBRDatasetSidecar.cc (checks F1-F11 and F13, codes
BBR024 and BBR025) and of the sidecar checks in BBRHFSSData (C1-C5, on both
polarizations), plus what the C++ does not do: the full mode lists of schema 1
(Bessel zeros through scipy), the CSV checksums, the far-field step, the agreement
of the two polarizations, and the builders the mock generators and the
legacy-sidecar script use. The schema (agreed with Blackbody-Simulations on
2026-10-05) is described in validation/README.md, Dataset sidecars. Errors are
ValueError whose message starts with the C++ code. invariant() and
same_invariant() hold what every frequency of one ID must share: the schema
major, the incidence grid and the physics, numbers to 1e-9 relative, over the
fields both sidecars carry; invariant_diff() names the blocks that differ.

Schema 2.x (OR-028, OR-031, OR-034, OR-036, OR-038, HF-040, HF-047, HF-049): the
checks run in the C++ order so the first error carries the same code. The 2.0
sidecar declares T = OutgoingPower / (IngoingPower * InterceptFactor), the per-key
divisor c_k from the declared incidence grid (cell_90 at theta = 90 deg), the
per-key t_flux, scale and far-field records, the quadrature and sampler, the
cross-term limit and the regime cutoffs; the waveguide table has 16 columns read by
position. check_csvs() adds the row checks of BBRHFSSData::LoadWaveguide (BBR012).
build_schema2_from_csvs() writes a 2.0 sidecar from finalized 2.0 tables.

Frames: canonical (p, l, g) = (propagation, long, gap), p x l = g, equal to the
Geant4 crack-local (x, y, z). The sampler implements one frame: the HFSS global
axes are X = -g, Y = +l, Z = +p, and the exit CS (Z, Y, -X) is (p, l, g).
"""
import hashlib
import json
import math
import os

import numpy as np
import pandas as pd

from . import hfss

C = 299792458.0                 # m/s, exact
X11_PRIME = 1.8411837813        # first zero of J1': the TE11 cutoff of a round guide
LIST_LIMIT_GHZ = 20000.0        # top of BBRsim's Planck band; list entries are <= this
TIE_REL = 1e-12                 # cutoffs this close are one tie: TE before TM, then by index
AXIS_TOL = 1e-9
MODE_REL = 1e-6
SECTION_REL = 1e-6              # section containment tolerance: the HFSS runner's; never stricter
POLARIZATION = "Ephi=0: E_theta=1; Ephi=1: E_phi=1"
TRANSMITTANCE = "outgoing_power_w / incoming_power_w, first row per key"
FF_COLUMNS = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta",
              "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"]
WG_COLUMNS = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower",
              "X", "Y", "Z", "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]
# The 2.0 list: headers in the stored (interleaved) order, InterceptFactor last. Read by position:
# Ex = (9, 10), Ey = (11, 12), Ez = (13, 14), InterceptFactor = 15 (BBRHFSSData's 2.x branch).
WG_COLUMNS_2 = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower",
                "X", "Y", "Z", "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor"]
# The fixed 2.0 strings (a differing value is BBR025), copied from the field list.
SCHEMA2_STRINGS = {
    "model.kind": "screen_free_abc",
    "transmittance.definition": "OutgoingPower / (IngoingPower * InterceptFactor)",
    "transmittance.numerator": "exit_flux_re_s_dot_n",
    "transmittance.incoming": "area_times_intensity_undivided",
    "transmittance.above_one": "per_photon_cap",
    "exit_field.sampler": "abs_e_t_squared_times_quadrature",
    "exit_field.scale.rule": "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w",
    "exit_field.column_order_check.rule": "reference_exporter_interleaved_v1",
    "cross_term.density": "e_t_quadrature",
    "cross_term.above_limit": "warn_once",
    "far_field.normalization": "pattern_only",
}
SCHEMA2_QUADRATURE = {"rectangle": "trapezoid", "disc": "uniform_disc"}
SCHEMA2_TANGENTIAL = ["Ex", "Ey"]
SCHEMA2_PLANE_TOL_M = 1e-12
SCHEMA2_CROSS_TERM = ("Re rho = sum_i w_i Re(E_t0,i . conj E_t1,i) / sqrt(sum_i w_i abs(E_t0,i)^2 * "
                      "sum_i w_i abs(E_t1,i)^2)")
KEY_TOL = 1e-9                  # per-key records match the declared keys (degrees)
DIVISOR_REL = 1e-6              # c_k, cell_90 against the closed form; the two T sources
STEP_REL = 1e-6
COS_FLOOR = 0.01
DARK_REL = 1e-12                # the producer's threshold: |flux| <= DARK_REL * (A.I) is a numerical zero
X01 = 2.4048255577              # first zero of J0: TM01 of a round guide
X21_PRIME = 3.0542369282        # first zero of J2': TE21 of a round guide
ETA0 = 1.25663706127e-6 * C     # ohm; the HFSS side's constants (recorded values only)
HFSS_AXES = {"x": [0, 0, -1], "y": [0, 1, 0], "z": [1, 0, 0]}        # in canonical components
EXIT_CS_GLOBAL = {"x": [0, 0, 1], "y": [0, 1, 0], "z": [-1, 0, 0]}   # in HFSS global components
IDENTITY = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
CANONICAL = "p,l,g; p x l = g; Geant4 crack-local (x,y,z) = (p,l,g)"   # the HFSS writer's exact text


def path_for(base, stem):
    return os.path.join(base, f"{stem}.dataset.json")


def load(base, stem):
    p = path_for(base, stem)
    try:
        with open(p, encoding="utf-8-sig") as fh:
            return json.load(fh)
    except FileNotFoundError:
        raise ValueError(f"BBR024: {p}: missing; every HFSS dataset needs this sidecar") from None
    except json.JSONDecodeError as e:
        raise ValueError(f"BBR024: {p}: not valid JSON ({e})") from None


def write(base, stem, sc):
    p = path_for(base, stem)
    with open(p, "w") as fh:
        json.dump(sc, fh, indent=2)
        fh.write("\n")
    return p


# --- modes -------------------------------------------------------------------

def mode_label(kind, i, j):
    return f"{kind}{i}{j}" if i < 10 and j < 10 else f"{kind}{i},{j}"


def _rect_f(a_m, b_m, m, n):
    return C / 2 * math.hypot(m / a_m, n / b_m) / 1e9


def _rect_count(a_m, b_m, limit_ghz):
    """TE (m, n) != (0, 0) plus TM (m, n >= 1) with cutoff <= limit."""
    total = 0
    for n in range(int(2 * b_m * limit_ghz * 1e9 / C) + 2):
        for m in range(int(2 * a_m * limit_ghz * 1e9 / C) + 2):
            if (m, n) == (0, 0) or _rect_f(a_m, b_m, m, n) > limit_ghz:
                continue
            total += 2 if (m >= 1 and n >= 1) else 1
    return total


def _rect_lowest(a_m, b_m):
    """(mode, cutoff_ghz) of the lowest mode of a closed rectangular guide: TE10 if a >= b, else TE01."""
    return ("TE10", _rect_f(a_m, b_m, 1, 0)) if a_m >= b_m else ("TE01", _rect_f(a_m, b_m, 0, 1))


def rect_modes(a_m, b_m, freq_ghz, limit_ghz=LIST_LIMIT_GHZ):
    """Closed ideal-PEC rectangular guide, a along l (index m), b along g (index n)."""
    mode, cut = _rect_lowest(a_m, b_m)
    onsets = [{"n": 0, "mode": "TE10", "cutoff_ghz": _rect_f(a_m, b_m, 1, 0)}]
    n = 1
    while _rect_f(a_m, b_m, 0, n) <= limit_ghz:
        onsets.append({"n": n, "mode": mode_label("TE", 0, n), "cutoff_ghz": _rect_f(a_m, b_m, 0, n)})
        n += 1
    return {"cutoff_ghz": cut, "mode": mode, "list_limit_ghz": limit_ghz,
            "polarization_filter_limit_ghz": _rect_f(a_m, b_m, 0, 1),
            "gap_family_onsets": onsets,
            "mode_count_below_limit": _rect_count(a_m, b_m, limit_ghz),
            "propagating_count": _rect_count(a_m, b_m, freq_ghz)}


def _order(entries):
    """Ascending cutoff; a run of cutoffs within TIE_REL of its first is sorted TE, TM, then index."""
    entries = sorted(entries)
    out, i = [], 0
    while i < len(entries):
        j = i + 1
        while j < len(entries) and entries[j][0] - entries[i][0] <= TIE_REL * entries[i][0]:
            j += 1
        out += sorted(entries[i:j], key=lambda e: (e[1] != "TE", e[2], e[3]))
        i = j
    return out


def _disc_entries(radius_m, limit_ghz):
    """(cutoff_ghz, kind, m, p) of every TE_mp and TM_mp of a closed PEC circular guide with its cutoff
    at or below limit_ghz, in list order."""
    from scipy.special import jn_zeros, jnp_zeros
    scale = C / (2 * math.pi * radius_m) / 1e9         # GHz per unit Bessel zero
    x_max = limit_ghz / scale
    entries, m = [], 0
    while True:
        added = False
        for kind, zeros in (("TE", jnp_zeros), ("TM", jn_zeros)):
            k = 8
            z = zeros(m, k)
            while z[-1] <= x_max:
                k *= 2
                z = zeros(m, k)
            for p, x in enumerate(z, start=1):
                if x <= x_max:
                    entries.append((float(x) * scale, kind, m, p))
                    added = True
        # Order 0 can be empty while TE11 lies below the limit (x'_11 < j_01 < x'_01); from order 1
        # on the first zero grows with the order, so the first empty order ends the search (CDX-006).
        if not added and m >= 1:
            break
        m += 1
    return _order(entries)


def disc_modes(radius_m, freq_ghz, limit_ghz=LIST_LIMIT_GHZ):
    """Closed ideal-PEC circular guide: TE_np = x'_np c/(2 pi R), TM_np = x_np c/(2 pi R). The list
    stops at limit_ghz; propagating_count counts to the file's own frequency, above the limit too
    (CDX-007), as the HFSS producer does."""
    cutoffs = [{"mode": mode_label(k, m, p), "cutoff_ghz": f, "degeneracy": 1 if m == 0 else 2}
               for f, k, m, p in _disc_entries(radius_m, limit_ghz)]
    lowest = X11_PRIME * C / (2 * math.pi * radius_m) / 1e9
    return {"cutoff_ghz": lowest, "mode": "TE11", "list_limit_ghz": limit_ghz,
            "polarization_filter_limit_ghz": lowest, "cutoffs": cutoffs,
            "propagating_count": len(_disc_entries(radius_m, freq_ghz))}


def modes_for(section, freq_ghz):
    shape = section["shape"]
    if shape == "rectangle":
        a, b = 2 * section["y_e_half_m"], 2 * section["z_e_half_m"]
        m = rect_modes(a, b, freq_ghz)
        m["basis"] = (f"closed ideal-PEC rectangular guide, a = {a!r} m along l, b = {b!r} m along g; "
                      "f_mn = (c/2) sqrt((m/a)^2 + (n/b)^2); no TEM mode")
        return m
    if shape == "disc":
        r = section["radius_m"]
        m = disc_modes(r, freq_ghz)
        m["basis"] = (f"closed ideal-PEC circular guide, R = {r!r} m; "
                      "TE_np = x'_np c/(2 pi R), TM_np = x_np c/(2 pi R)")
        return m
    raise ValueError(f"BBR025: cross_section shape {shape!r} has no mode model")


def lowest_mode(section):
    """(mode, cutoff_ghz, polarization_filter_limit_ghz): what the C++ F13 re-derives."""
    if section["shape"] == "rectangle":
        a, b = 2 * section["y_e_half_m"], 2 * section["z_e_half_m"]
        mode, cut = _rect_lowest(a, b)
        return mode, cut, _rect_f(a, b, 0, 1)
    if section["shape"] == "disc":
        cut = X11_PRIME * C / (2 * math.pi * section["radius_m"]) / 1e9
        return "TE11", cut, cut
    raise ValueError(f"BBR025: cross_section shape {section['shape']!r} has no mode model")


# --- schema 2: the divisor, the regime cutoffs --------------------------------

def schema_major(sc):
    """1 for schema_version "1" or "1.x", 2 for "2" or "2.x"; anything else is BBR024 (the version gate)."""
    v = _req(sc, "schema_version", "", "str")
    if v == "1" or v.startswith("1."):
        return 1
    if v == "2" or v.startswith("2."):
        return 2
    raise ValueError(f"BBR024: schema_version {v}: BBRsim reads 1.x and 2.x")


def _delta(theta_list):
    """(smallest declared theta above 90 deg) - 90, or -1 without a 90 deg entry or anything above it."""
    has90, above = False, None
    for t in theta_list:
        if abs(t - 90.0) <= KEY_TOL:
            has90 = True
        elif t > 90.0:
            above = t if above is None else min(above, t)
    return above - 90.0 if has90 and above is not None else -1.0


def cell_90(theta_list):
    """HF-049 / OR-036: 1 / <1 / max(cos a, 0.01)> with a uniform on [90 - D/2, 90] deg, D = (smallest
    declared theta above 90) - 90. None for a declared list without 90 deg or without anything above it."""
    delta = _delta(theta_list)
    if delta < 0.0:
        return None
    a0 = math.pi / 2 - delta / 2 * math.pi / 180.0
    af = math.acos(COS_FLOOR)
    length = math.pi / 2 - a0

    def prim(a):                                  # integral of sec
        return math.log(1.0 / math.cos(a) + math.tan(a))

    inv = (prim(af) - prim(a0)) + (math.pi / 2 - af) / COS_FLOOR if af > a0 else (math.pi / 2 - a0) / COS_FLOOR
    return length / inv


def intercept_factor(theta_deg, theta_list):
    """c_k: |cos theta|, or cell_90 at theta = 90 deg (None if the list gives no cell)."""
    if abs(theta_deg - 90.0) <= KEY_TOL:
        return cell_90(theta_list)
    return abs(math.cos(theta_deg * math.pi / 180.0))


def regime_cutoffs(section):
    """The cutoffs a photon can cross to reach a table (OR-031), ascending: a rectangle's TE10 and gap
    family TE0n up to 20 THz inclusive; a disc's TE11, TM01 and TE21 (closed-form zeros, no scipy).
    A list of {"mode", "cutoff_ghz"}."""
    if section["shape"] == "disc":
        k = C / (2 * math.pi * section["radius_m"]) / 1e9
        return [{"mode": "TE11", "cutoff_ghz": X11_PRIME * k}, {"mode": "TM01", "cutoff_ghz": X01 * k},
                {"mode": "TE21", "cutoff_ghz": X21_PRIME * k}]
    if section["shape"] != "rectangle":
        raise ValueError(f"BBR025: cross_section shape {section['shape']!r} has no mode model")
    a, b = 2 * section["y_e_half_m"], 2 * section["z_e_half_m"]
    out = [{"mode": "TE10", "cutoff_ghz": C / (2 * a) / 1e9}]
    n = 1
    while n * C / (2 * b) / 1e9 <= LIST_LIMIT_GHZ + 1e-9:
        out.append({"mode": f"TE0{n}" if n < 10 else f"TE0,{n}", "cutoff_ghz": n * C / (2 * b) / 1e9})
        n += 1
    return sorted(out, key=lambda e: e["cutoff_ghz"])          # stable, as std::stable_sort


# --- CSV description, checksums, builder --------------------------------------

def _range(values):
    v = np.unique(np.asarray(values, float))
    return {"min": float(v[0]), "max": float(v[-1]), "count": int(len(v))}


FIELD_COLUMNS = WG_COLUMNS[9:15]   # the six exit-field components (real and imaginary E)


def _read_csv(base, stem, ephi, name, columns):
    """One CSV as a DataFrame; a missing file is BBR001 (far field) or BBR002 (waveguide), as in
    BBRHFSSData, and a header other than the positional order BBRHFSSData reads is BBR013 (C1).
    A data row with the wrong number of fields is BBR013, the loader's code for a malformed row."""
    p = os.path.join(base, f"{stem}_Ephi={ephi}", name)
    try:
        df = pd.read_csv(p)
    except FileNotFoundError:
        raise ValueError(f"{'BBR001' if name == 'far_field.csv' else 'BBR002'}: {p}: cannot open") from None
    except pd.errors.EmptyDataError:
        raise ValueError(f"BBR013: {p}: empty, no header") from None
    except pd.errors.ParserError as e:
        raise ValueError(f"BBR013: {p}: a row has the wrong number of fields ({str(e).strip()})") from None
    if list(df.columns) != columns:
        raise ValueError(f"BBR013: {p}: header {list(df.columns)} differs from the columns "
                         f"BBRHFSSData reads {columns}")
    # pandas takes the first field as an index when every data row has one field more than the
    # header, and fills a short row with NaN; both are rows with the wrong number of fields
    if not isinstance(df.index, pd.RangeIndex):
        raise ValueError(f"BBR013: {p}: data rows have more fields than the header")
    na = df.isna().to_numpy()
    if na.any():
        r, c = (int(x) for x in np.argwhere(na)[0])
        raise ValueError(f"BBR013: {p}: a row lacks fields or holds an empty field or NaN "
                         f"(line {r + 2}, column {columns[c]})")
    # every field BBRHFSSData parses with Num, all but Freq and Ephi, is a finite number (CDX-001)
    numeric = [c for c in columns if c not in ("Freq", "Ephi")]
    bad = ~np.isfinite(df[numeric].apply(pd.to_numeric, errors="coerce").to_numpy(float))
    if bad.any():
        r, c = (int(x) for x in np.argwhere(bad)[0])
        raise ValueError(f"BBR013: {p}: line {r + 2}, column {numeric[c]}: '{df[numeric[c]].iloc[r]}' is not "
                         f"a finite number")
    return df


def _rows_by_key(df, cols):
    """{key: array of the cols of each incidence key's rows, in file order}, keys rounded to 0.01 deg
    as BBRHFSSData rounds them."""
    out = {}
    for (p, t), g in df.groupby(["IWavePhi", "IWaveTheta"], sort=False):
        k = (round(float(p), 2), round(float(t), 2))
        v = g[list(cols)].to_numpy(float)
        out[k] = np.vstack([out[k], v]) if k in out else v
    return out


def _describe_one(base, stem, ephi, v2=False):
    """What the two CSVs of one polarization contain: header, incidence keys, grids, per-key counts.
    v2: the waveguide table is the 2.0 one (16 columns); the row arrays the 2.x row checks need are kept."""
    ff = _read_csv(base, stem, ephi, "far_field.csv", FF_COLUMNS)
    wg = _read_csv(base, stem, ephi, "waveguide.csv", WG_COLUMNS_2 if v2 else WG_COLUMNS)
    field_cols = WG_COLUMNS_2[9:15] if v2 else FIELD_COLUMNS

    def per_key(df, what):
        n = df.groupby(["IWavePhi", "IWaveTheta"]).size().unique()
        if len(n) != 1:
            raise ValueError(f"BBR012: {stem}_Ephi={ephi}: {what} rows per incidence key differ ({sorted(n)})")
        return int(n[0])

    th, ph = _range(ff.Theta), _range(ff.Phi)
    for r in (th, ph):
        r["step"] = (r["max"] - r["min"]) / (r["count"] - 1) if r["count"] > 1 else 0.0
    rows_key = [(round(float(p), 2), round(float(t), 2)) for p, t in
                wg[["IWavePhi", "IWaveTheta"]].itertuples(index=False)]
    yz = wg[["Y", "Z"]].to_numpy(float)
    return {
        "ephi": ephi,
        "incident_phi_deg": sorted(float(v) for v in wg.IWavePhi.unique()),
        "incident_theta_deg": sorted(float(v) for v in wg.IWaveTheta.unique()),
        "far_field": {"theta_deg": th, "phi_deg": ph, "points_per_key": per_key(ff, "far-field"),
                      "columns": list(ff.columns)},
        "exit_field": {"grid": {"y_e": _range(wg.Y), "z_e": _range(wg.Z)},
                       "points_per_key_retained": per_key(wg, "exit-point"), "columns": list(wg.columns)},
        "incoming_power_w": float(wg.IngoingPower.iloc[0]),
        "x_max_abs_m": float(np.abs(wg.X.to_numpy(float)).max()),
        "yz_m": yz,
        # (key, Y, Z) of every exit row whose six field components are exactly zero
        "zero_points": {(k, float(y), float(z)) for k, (y, z), zero in
                        zip(rows_key, yz, (wg[field_cols].to_numpy(float) == 0).all(axis=1)) if zero},
        "rows_key": rows_key,
        "keys": set(rows_key),
        "ff_keys": {(round(float(p), 2), round(float(t), 2)) for p, t in
                    ff[["IWavePhi", "IWaveTheta"]].itertuples(index=False)},
        "freq_tokens": {"far_field.csv": sorted(ff.Freq.astype(str).unique()),
                        "waveguide.csv": sorted(wg.Freq.astype(str).unique())},
        "ff_rows": _rows_by_key(ff, ("Phi", "Theta")),
        "exit_rows": _rows_by_key(wg, ("X", "Y", "Z")),
        # 2.x, in file order: the raw (IWavePhi, IWaveTheta) of each row, OutgoingPower, IngoingPower and
        # InterceptFactor, for the row checks of BBRHFSSData::LoadWaveguide
        "rows2": ((wg[["IWavePhi", "IWaveTheta"]].to_numpy(float), wg.OutgoingPower.to_numpy(float),
                   wg.IngoingPower.to_numpy(float), wg.InterceptFactor.to_numpy(float)) if v2 else None),
    }


def describe_csvs(base, stem, v2=False):
    """What the CSVs of <stem> contain. The top-level fields describe Ephi=0 (what the builder
    records); "polarizations" holds the description of each of Ephi=0 and Ephi=1. v2: 2.0 tables."""
    d0, d1 = _describe_one(base, stem, 0, v2), _describe_one(base, stem, 1, v2)
    return {**d0, "polarizations": (d0, d1)}


def file_entries(base, stem):
    """{"<stem>_Ephi=N/<file>": {sha256, bytes, rows}} for the four CSVs. rows is the number of
    non-blank lines after the header (a blank line holds only whitespace), which is pandas'
    len(read_csv(...)) for these files, with or without a trailing newline."""
    out = {}
    for e in (0, 1):
        for name in ("far_field.csv", "waveguide.csv"):
            rel = f"{stem}_Ephi={e}/{name}"
            p = os.path.join(base, rel)
            with open(p, "rb") as fh:
                data = fh.read()
            lines = sum(1 for line in data.splitlines() if line.strip())
            out[rel] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data),
                        "rows": max(lines - 1, 0)}
    return out


def build_from_csvs(base, dataset_id, frequency_label, *, section, extent_mm, provenance,
                    exit_origin_mm, plane_wave_origin_mm, bounding_box_mm, pose_rule="canonical-z",
                    resolution_mm=None, outside_points="none", rim_points="included",
                    bounds_method="vertices", rotational=False, walls="PEC"):
    """A schema-1.0 sidecar for <dataset_id>_<frequency_label>, read off its CSVs."""
    stem = f"{dataset_id}_{frequency_label}"
    freq = hfss.parse_frequency_GHz(frequency_label)
    d = describe_csvs(base, stem)
    return {
        "schema_version": "1.0",
        "dataset_id": dataset_id,
        "frequency_ghz": freq,
        "frequency_label": frequency_label,
        "provenance": provenance,
        "frames": {
            "canonical": CANONICAL, "pose_rule": pose_rule,
            "hfss_global_axes_in_canonical": HFSS_AXES,
            "entrance_face": {"axis": "z", "side": "min"}, "exit_face": {"axis": "z", "side": "max"},
            "entrance_outward_normal_global": [0, 0, -1], "exit_outward_normal_global": [0, 0, 1],
            "exit_cs": {"name": "outgoing_cs", "origin": "exit_face_center",
                        "origin_mm_global": exit_origin_mm, **EXIT_CS_GLOBAL},
            "exit_cs_axes_in_canonical": IDENTITY,
        },
        "excitation": {
            "coordinate_system": "global", "incidence_convention": "arrival_direction",
            "normal_entry_theta_deg": 180.0,
            "incident_phi_deg": d["incident_phi_deg"], "incident_theta_deg": d["incident_theta_deg"],
            "ei_v_per_m": 1.0, "incoming_power_w": d["incoming_power_w"],
            "polarization_convention": POLARIZATION,
            "plane_wave_origin": "entrance_face_center", "origin_mm_global": plane_wave_origin_mm,
        },
        "far_field": {"coordinate_system": "exit_cs", "definition": "Theta-Phi",
                      "component_basis": "spherical_in_exit_cs", "radiation_surface": "exit_face",
                      **d["far_field"]},
        "exit_field": {
            "coordinate_system": "exit_cs", "points_in_si": True, "field_in_ref_cs": False,
            "field_components_frame": "hfss_global", "plane": "x_e=0", "resolution_mm": resolution_mm,
            "grid": d["exit_field"]["grid"], "cross_section": section, "outside_points": outside_points,
            "bounds_method": bounds_method, "rim_points": rim_points,
            "points_per_key_retained": d["exit_field"]["points_per_key_retained"],
            "columns": d["exit_field"]["columns"],
        },
        "transmittance": {"definition": TRANSMITTANCE, "outgoing_power": "integral |Re S| over exit face",
                          "incoming_includes_cos_theta": False},
        "symmetry": {"mirror_l": True, "mirror_g": True, "end_to_end": True, "rotational": rotational},
        "boundaries": {"entrance": "radiation", "exit": "radiation", "walls": walls},
        "geometry": {"shape": "box" if section["shape"] == "rectangle" else "cylinder",
                     "extent_mm": extent_mm, "bounding_box_mm": bounding_box_mm},
        "modes": modes_for(section, freq),
        "files": file_entries(base, stem),
    }


INVARIANT_REL = 1e-9


def invariant(sc):
    """The frequency-independent physics of a sidecar (BBRDatasetSidecar::invariant).

    Every frequency of one dataset must share it, whichever writer or pose produced
    each file: the schema major and the incidence grid (a tree holds one of each per
    dataset ID, G4-032), the frame mapping, symmetry, boundaries, geometry shape and
    extents, the modes block without basis (recorded only) and propagating_count (per
    frequency), and the exit cross-section. Descriptive and pose-dependent fields
    (canonical text, pose_rule, face selectors, origins, bounding box) are left out.
    A block the sidecar lacks is left out of the comparison, as _close() compares the keys both carry.
    """
    fr, geo = sc["frames"], sc["geometry"]
    inv = {}
    if "schema_version" in sc:
        inv["schema"] = schema_major(sc)
    if "excitation" in sc:
        inv["excitation"] = {k: sc["excitation"][k] for k in ("incident_phi_deg", "incident_theta_deg")
                             if k in sc["excitation"]}
    inv.update({
        "frames": {k: fr[k] for k in ("hfss_global_axes_in_canonical", "exit_cs_axes_in_canonical")},
        "symmetry": sc["symmetry"], "boundaries": sc["boundaries"],
        "geometry": {k: geo[k] for k in ("shape", "extent_mm") if k in geo},
        "modes": {k: v for k, v in sc["modes"].items() if k not in ("propagating_count", "basis")},
        "cross_section": sc["exit_field"]["cross_section"],
    })
    return inv


def _close(x, y):
    """Objects over the keys both carry, lists element-wise, numbers to INVARIANT_REL, the rest exactly."""
    if isinstance(x, dict) and isinstance(y, dict):
        return all(_close(x[k], y[k]) for k in x.keys() & y.keys())
    if isinstance(x, list) and isinstance(y, list):
        return len(x) == len(y) and all(_close(a, b) for a, b in zip(x, y))
    num = (int, float)
    if isinstance(x, num) and isinstance(y, num) and not isinstance(x, bool) and not isinstance(y, bool):
        return abs(x - y) <= INVARIANT_REL * max(abs(x), abs(y))
    return x == y


def same_invariant(a, b):
    """True when two sidecars of one dataset carry the same frequency-independent physics."""
    return _close(invariant(a), invariant(b))


def invariant_diff(a, b):
    """The top-level invariant() blocks on which two sidecars differ, in invariant() order;
    empty exactly when same_invariant(a, b) holds."""
    ia, ib = invariant(a), invariant(b)
    return [k for k in ia if k in ib and not _close(ia[k], ib[k])]


# --- checks --------------------------------------------------------------------

def _req(d, key, where, kind=None):
    if not isinstance(d, dict) or key not in d:
        raise ValueError(f"BBR024: missing field {where}{key}")
    v = d[key]
    num = isinstance(v, (int, float)) and not isinstance(v, bool)
    ok = {None: True, "str": isinstance(v, str), "bool": isinstance(v, bool),
          "num": num and math.isfinite(v), "int": isinstance(v, int) and not isinstance(v, bool),
          "dict": isinstance(v, dict), "list": isinstance(v, list)}[kind]
    if not ok:
        raise ValueError(f"BBR024: {where}{key} has the wrong type ({type(v).__name__})")
    return v


def _expect(got, want, field):
    if got != want:
        raise ValueError(f"BBR025: {field} is {got!r}; BBRsim implements only {want!r}")


def _vec(d, key, where):
    v = _req(d, key, where, "list")
    if len(v) != 3 or not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
        raise ValueError(f"BBR024: {where}{key} must be a 3-vector of numbers")
    return np.asarray(v, float)


def _near(a, b):
    return bool(np.all(np.abs(np.asarray(a, float) - np.asarray(b, float)) <= AXIS_TOL))


def _axis_range(d, key, where):
    r = _req(d, key, where, "dict")
    w = f"{where}{key}."
    lo, hi, n = _req(r, "min", w, "num"), _req(r, "max", w, "num"), _req(r, "count", w, "int")
    if n < 1 or hi < lo:
        raise ValueError(f"BBR024: {w} needs min <= max and count >= 1")
    return lo, hi, n


def _rep(v):
    """A number as BBRsim's messages print it: the shortest text that reads back, no trailing .0."""
    v = float(v)
    return str(int(v)) if v == int(v) and abs(v) < 1e15 else repr(v)


def _g12(v):
    return f"{float(v):.12g}"


def _ks(phi, theta):
    return f"(phi {_g12(phi)}, theta {_g12(theta)})"


def _relnear(a, b, rel):
    return abs(a - b) <= rel * abs(b)


def _nums(d, key, where):
    v = _req(d, key, where, "list")
    if not v or not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in v):
        raise ValueError(f"BBR024: {where}{key} must be a non-empty array of numbers")
    return [float(x) for x in v]


def _match_keys(parent, key, where, phis, thetas):
    """The records of one per-key list, phi-major then theta, each matched to a declared key within
    KEY_TOL; a duplicate, an unmatched record or a gap is BBR024 (MatchKeys)."""
    lst = _req(parent, key, where)
    w = where + key
    if not isinstance(lst, list):
        raise ValueError(f"BBR024: {w} must be a list of per-key records")
    out = [None] * (len(phis) * len(thetas))
    for i, rec in enumerate(lst):
        wi = f"{w}[{i}]."
        if not isinstance(rec, dict):
            raise ValueError(f"BBR024: {w}[{i}] must be an object")
        phi, theta = _req(rec, "phi_deg", wi, "num"), _req(rec, "theta_deg", wi, "num")
        ip = next((a for a, x in enumerate(phis) if abs(phi - x) <= KEY_TOL), None)
        it = next((b for b, x in enumerate(thetas) if abs(theta - x) <= KEY_TOL), None)
        if ip is None or it is None:
            raise ValueError(f"BBR024: {wi}key {_ks(phi, theta)} does not match a declared key of "
                             "excitation.incident_phi_deg x incident_theta_deg")
        slot = ip * len(thetas) + it
        if out[slot] is not None:
            raise ValueError(f"BBR024: {w}: duplicate record for key {_ks(phis[ip], thetas[it])}")
        out[slot] = rec
    for ip, phi in enumerate(phis):
        for it, theta in enumerate(thetas):
            if out[ip * len(thetas) + it] is None:
                raise ValueError(f"BBR024: {w} has no record for declared key {_ks(phi, theta)}")
    return out


def _pols(rec, where):
    """The two polarizations of a per-key record: exactly "0" and "1" (Polarizations)."""
    pols = _req(rec, "polarizations", where)
    if not isinstance(pols, dict) or len(pols) != 2:
        raise ValueError(f'BBR024: {where}polarizations must be an object with exactly "0" and "1"')
    return _req(pols, "0", where + "polarizations."), _req(pols, "1", where + "polarizations.")


def _schema2_keys(sc, phis, thetas):
    """The per-key records of a 2.x sidecar (BBRKeyRecord), in declared order, after every per-key
    check of ParseSchema2: c_key against the closed form (BBR025) and per_key (BBR024), the scale
    weights, the dark flags against t_flux (OR-038), the far-field ratios."""
    tr = _req(sc, "transmittance", "", "dict")
    c_recs = _match_keys(tr, "c_key", "transmittance.", phis, thetas)
    k_recs = _match_keys(tr, "per_key", "transmittance.", phis, thetas)
    s_recs = _match_keys(_req(_req(sc, "exit_field", "", "dict"), "scale", "exit_field.", "dict"), "per_key",
                         "exit_field.scale.", phis, thetas)
    f_recs = _match_keys(_req(sc, "far_field", "", "dict"), "ff_over_flux", "far_field.", phis, thetas)
    out = []
    for ip, phi in enumerate(phis):
        for it, theta in enumerate(thetas):
            k = ip * len(thetas) + it
            ks = _ks(phi, theta)
            ck = _req(c_recs[k], "c_k", f"transmittance.c_key {ks} ", "num")
            want = intercept_factor(theta, thetas)
            if not _relnear(ck, want, DIVISOR_REL):
                raise ValueError(f"BBR025: transmittance.c_key {ks} c_k is {_rep(ck)}, but the declared grid gives "
                                 f"{_rep(want)} (|cos theta|, or cell_90 at 90 deg; 1e-6 relative)")
            pk = k_recs[k]
            wp = f"transmittance.per_key {ks} "
            pk_c = _req(pk, "c_k", wp, "num")
            if not _relnear(pk_c, ck, DIVISOR_REL):
                raise ValueError(f"BBR024: {wp}c_k {_rep(pk_c)} disagrees with transmittance.c_key {_rep(ck)}")
            sk = s_recs[k]
            ws = f"exit_field.scale.per_key {ks} "
            w0, w1 = _req(sk, "w0", ws, "num"), _req(sk, "w1", ws, "num")
            if not (w0 > 0.0 and w1 > 0.0):
                raise ValueError(f"BBR024: {ws}w0 and w1 must be > 0 (w0 {_rep(w0)}, w1 {_rep(w1)})")
            dark = (_req(sk, "dark0", ws, "bool"), _req(sk, "dark1", ws, "bool"))
            tf = []
            for p, pol in enumerate(_pols(pk, wp)):
                wpp = f"{wp}polarizations.{p}."
                t = _req(pol, "t_flux", wpp, "num")
                _req(pol, "outgoing_power_abs_re_s_w", wpp, "num")
                if dark[p] and t != 0.0:
                    raise ValueError(f"BBR024: key {ks} polarization {p} is dark (exit_field.scale.per_key dark{p}), "
                                     f"but transmittance.per_key t_flux is {_g12(t)}; a dark basis needs t_flux "
                                     "exactly 0 (OR-038)")
                tf.append(t)
            for p, v in enumerate(_pols(f_recs[k], f"far_field.ff_over_flux {ks} ")):
                wf = f"far_field.ff_over_flux {ks} polarizations.{p}"
                if v is None:
                    if not dark[p]:
                        raise ValueError(f"BBR024: {wf} is null, but polarization {p} of key {ks} is not dark")
                elif isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                    raise ValueError(f"BBR024: {wf} must be a finite number, or null for a dark basis")
            out.append({"phi": phi, "theta": theta, "c_k": ck, "w0": w0, "w1": w1, "dark0": dark[0], "dark1": dark[1],
                        "t_flux0": tf[0], "t_flux1": tf[1]})
    return out


def key_records(sc):
    """The 2.x per-key records of a sidecar that has passed check(): a dict per declared key in declared
    order (phi-major, then theta) with phi, theta, c_k, w0, w1, dark0, dark1, t_flux0, t_flux1."""
    exc = sc["excitation"]
    return _schema2_keys(sc, [float(x) for x in exc["incident_phi_deg"]], [float(x) for x in exc["incident_theta_deg"]])


def _schema2(sc, section):
    """ParseSchema2: what a 2.x sidecar declares in place of the 1.x transmittance block."""
    exc = sc["excitation"]
    phis, thetas = [float(x) for x in exc["incident_phi_deg"]], [float(x) for x in exc["incident_theta_deg"]]
    S = SCHEMA2_STRINGS
    _expect(_req(_req(sc, "model", "", "dict"), "kind", "model.", "str"), S["model.kind"], "model.kind")
    _req(sc, "conventions", "", "dict")

    tr = _req(sc, "transmittance", "", "dict")                                  # OR-034, HF-049
    want90 = cell_90(thetas)
    if want90 is None:
        raise ValueError("BBR025: excitation.incident_theta_deg has no theta = 90 deg entry with a declared theta "
                         "above it, so transmittance.cell_90 cannot be derived from the declared grid (HF-049); "
                         "BBRsim refuses it")
    delta = _delta(thetas)
    for key in ("definition", "numerator", "incoming", "above_one"):
        _expect(_req(tr, key, "transmittance.", "str"), S[f"transmittance.{key}"], f"transmittance.{key}")
    cell = _req(tr, "cell_90", "transmittance.", "num")
    if not _relnear(cell, want90, DIVISOR_REL):
        raise ValueError(f"BBR025: transmittance.cell_90 is {_rep(cell)}, but the declared grid's outermost band "
                         f"gives {_rep(want90)} (1e-6 relative; HF-049)")
    wa = "transmittance.cell_90_assumes."
    assumes = _req(tr, "cell_90_assumes", "transmittance.", "dict")
    _expect(_req(assumes, "emitter", wa, "str"), "uniform_theta", wa + "emitter")   # BBRsim's one emitter
    floor = _req(assumes, "cos_floor", wa, "num")
    if not _relnear(floor, COS_FLOOR, 1e-9):
        raise ValueError(f"BBR025: {wa}cos_floor is {_rep(floor)}; BBRsim implements only {_rep(COS_FLOOR)}")
    band = _nums(assumes, "band_deg", wa)
    if len(band) != 2:
        raise ValueError(f"BBR024: {wa}band_deg must be [lo, 90]")
    if abs(band[0] - (90.0 - delta / 2)) > KEY_TOL or abs(band[1] - 90.0) > KEY_TOL:
        raise ValueError(f"BBR025: {wa}band_deg is [{_rep(band[0])}, {_rep(band[1])}]; the declared grid gives "
                         f"[{_rep(90.0 - delta / 2)}, 90]")

    _schema2_keys(sc, phis, thetas)

    xf = sc["exit_field"]                                                       # sampler, quadrature, fitted lattice
    rectangle = section["shape"] == "rectangle"
    _expect(_req(xf, "quadrature", "exit_field.", "str"), SCHEMA2_QUADRATURE[section["shape"]], "exit_field.quadrature")
    _expect(_req(xf, "sampler", "exit_field.", "str"), S["exit_field.sampler"], "exit_field.sampler")
    if _req(xf, "tangential_components", "exit_field.", "list") != SCHEMA2_TANGENTIAL:
        raise ValueError('BBR025: exit_field.tangential_components must be ["Ex", "Ey"]; the sampler reads the '
                         "exit-face tangentials only")
    _expect(_req(_req(xf, "scale", "exit_field.", "dict"), "rule", "exit_field.scale.", "str"),
            S["exit_field.scale.rule"], "exit_field.scale.rule")
    co = _req(xf, "column_order_check", "exit_field.", "dict")
    wc = "exit_field.column_order_check."
    _expect(_req(co, "rule", wc, "str"), S["exit_field.column_order_check.rule"], wc + "rule")
    _req(co, "method", wc, "str")
    if not _req(_req(co, "reference", wc, "dict"), "passed", wc + "reference.", "bool"):
        raise ValueError(f"BBR025: {wc}reference.passed is false: the exporter's column order was not validated")
    tol = _req(xf, "plane_tolerance_m", "exit_field.", "num")
    if not _relnear(tol, SCHEMA2_PLANE_TOL_M, 1e-9):
        raise ValueError(f"BBR025: exit_field.plane_tolerance_m is {_rep(tol)}; BBRsim implements only 1e-12")
    for axis in ("y_e", "z_e"):
        wg = f"exit_field.grid.{axis}"
        ga = xf["grid"][axis]
        lo, hi, n = ga["min"], ga["max"], ga["count"]
        if rectangle and n < 2:
            raise ValueError(f"BBR024: {wg} needs count >= 2 on a rectangle: neither the lattice step nor the "
                             "trapezoid weights are defined on one node")
        step = _req(ga, "step", wg + ".", "num")
        if not step > 0.0:
            raise ValueError(f"BBR024: {wg}.step must be > 0")
        if n >= 2 and not _relnear(step, (hi - lo) / (n - 1), STEP_REL):
            raise ValueError(f"BBR024: {wg}.step is {_rep(step)}, but (max - min)/(count - 1) is {_rep((hi - lo) / (n - 1))}")

    ct = _req(sc, "cross_term", "", "dict")                                     # cross term and far field
    _expect(_req(ct, "density", "cross_term.", "str"), S["cross_term.density"], "cross_term.density")
    _req(ct, "definition", "cross_term.", "str")
    _req(ct, "valid_below_ghz", "cross_term.", "num")
    _expect(_req(ct, "above_limit", "cross_term.", "str"), S["cross_term.above_limit"], "cross_term.above_limit")
    _expect(_req(sc["far_field"], "normalization", "far_field.", "str"), S["far_field.normalization"],
            "far_field.normalization")


def _check_regime_cutoffs(sc, md, section, filter_limit_ghz):
    """modes.regime_cutoffs_ghz re-derived entry by entry (BBR025, 1e-6 relative), then cross_term's limit
    against modes.polarization_filter_limit_ghz."""
    lst = _req(md, "regime_cutoffs_ghz", "modes.")
    if not isinstance(lst, list):
        raise ValueError("BBR024: modes.regime_cutoffs_ghz must be an array of {mode, cutoff_ghz}")
    got = []
    for e in lst:
        if not isinstance(e, dict):
            raise ValueError("BBR024: modes.regime_cutoffs_ghz entries must be objects {mode, cutoff_ghz}")
        got.append((_req(e, "mode", "modes.regime_cutoffs_ghz[].", "str"),
                    _req(e, "cutoff_ghz", "modes.regime_cutoffs_ghz[].", "num")))
    want = [(e["mode"], e["cutoff_ghz"]) for e in regime_cutoffs(section)]
    for i in range(max(len(got), len(want))):
        at = f"modes.regime_cutoffs_ghz[{i}]"
        if i >= len(got):
            raise ValueError(f"BBR025: {at} is missing: the section gives {want[i][0]} at {_rep(want[i][1])} GHz "
                             f"({len(want)} entries, the sidecar declares {len(got)})")
        if i >= len(want):
            raise ValueError(f"BBR025: {at} is {got[i][0]} at {_rep(got[i][1])} GHz, but the section gives only "
                             f"{len(want)} entries")
        if got[i][0] != want[i][0] or not _relnear(got[i][1], want[i][1], MODE_REL):
            raise ValueError(f"BBR025: {at} is {got[i][0]} at {_rep(got[i][1])} GHz, but the section gives "
                             f"{want[i][0]} at {_rep(want[i][1])} GHz (1e-6 relative)")
    valid = sc["cross_term"]["valid_below_ghz"]
    if not _relnear(valid, filter_limit_ghz, MODE_REL):
        raise ValueError(f"BBR025: cross_term.valid_below_ghz is {_rep(valid)}, but modes.polarization_filter_limit_ghz "
                         f"is {_rep(filter_limit_ghz)} GHz; the cross term is valid below that limit (1e-6 relative)")


def check(sc, dataset_id, stem, dir_freq_ghz):
    """F1-F11 and F13, as BBRDatasetSidecar::Parse runs them."""
    if not isinstance(sc, dict):
        raise ValueError("BBR024: the top level must be a JSON object")
    major = schema_major(sc)                                                       # F1
    if _req(sc, "dataset_id", "", "str") != dataset_id:                           # F2
        raise ValueError(f"BBR024: dataset_id {sc['dataset_id']!r}, the directory says {dataset_id!r}")
    prefix = dataset_id + "_"                                                      # F3
    if not stem.startswith(prefix):
        raise ValueError(f"BBR024: directory stem {stem} does not start with {prefix}")
    label = _req(sc, "frequency_label", "", "str")
    if label != stem[len(prefix):]:
        raise ValueError(f"BBR024: frequency_label {label!r}, the directory says {stem[len(prefix):]!r}")
    f, from_label = _req(sc, "frequency_ghz", "", "num"), hfss.parse_frequency_GHz(label)
    if not dir_freq_ghz > 0:
        raise ValueError(f"BBR024: directory frequency {dir_freq_ghz} GHz is not positive")
    if from_label is None:
        raise ValueError(f"BBR024: frequency_label {label!r} does not parse as a frequency")
    if abs(from_label - dir_freq_ghz) > 1e-3 * dir_freq_ghz:
        raise ValueError(f"BBR024: frequency_label {label!r} ({from_label} GHz) disagrees with the "
                         f"directory frequency {dir_freq_ghz} GHz")
    if abs(f - dir_freq_ghz) > 1e-3 * dir_freq_ghz:
        raise ValueError(f"BBR024: frequency_ghz {f} disagrees with the directory frequency {dir_freq_ghz} GHz")

    fr = _req(sc, "frames", "", "dict")                                            # F4
    ax = _req(fr, "hfss_global_axes_in_canonical", "frames.", "dict")
    hx, hy, hz = (_vec(ax, a, "frames.hfss_global_axes_in_canonical.") for a in "xyz")
    M = np.stack([hx, hy, hz], axis=1)
    if not (np.allclose(M.T @ M, np.eye(3), rtol=0, atol=AXIS_TOL) and _near(np.cross(hx, hy), hz)):
        raise ValueError("BBR025: frames.hfss_global_axes_in_canonical is not a right-handed orthonormal basis")
    ecs = _req(fr, "exit_cs", "frames.", "dict")
    ex, ey, ez = (_vec(ecs, a, "frames.exit_cs.") for a in "xyz")
    if not _near(np.cross(ex, ey), ez):
        raise ValueError("BBR025: frames.exit_cs.z is not x cross y")
    ecc = _req(fr, "exit_cs_axes_in_canonical", "frames.", "dict")
    cx, cy, cz = (_vec(ecc, a, "frames.exit_cs_axes_in_canonical.") for a in "xyz")
    if not (_near(M @ ex, cx) and _near(M @ ey, cy) and _near(M @ ez, cz)):
        raise ValueError("BBR025: frames.exit_cs_axes_in_canonical disagrees with exit_cs mapped "
                         "through hfss_global_axes_in_canonical")
    if not (_near(M @ _vec(fr, "exit_outward_normal_global", "frames."), [1, 0, 0])
            and _near(M @ _vec(fr, "entrance_outward_normal_global", "frames."), [-1, 0, 0])):
        raise ValueError("BBR025: the exit face's outward normal must map to +p, the entrance face's to -p")
    if not all(_near(v, w) for v, w in ((hx, HFSS_AXES["x"]), (hy, HFSS_AXES["y"]), (hz, HFSS_AXES["z"]),
                                         (cx, IDENTITY["x"]), (cy, IDENTITY["y"]), (cz, IDENTITY["z"]))):
        raise ValueError("BBR025: the frames differ from the one the sampler implements "
                         "(HFSS X = -g, Y = +l, Z = +p; exit CS = (p, l, g))")

    exc = _req(sc, "excitation", "", "dict")                                       # F5
    _expect(_req(exc, "coordinate_system", "excitation.", "str"), "global", "excitation.coordinate_system")
    _expect(_req(exc, "incidence_convention", "excitation.", "str"), "arrival_direction",
            "excitation.incidence_convention")
    if abs(_req(exc, "normal_entry_theta_deg", "excitation.", "num") - 180.0) > 1e-9:
        raise ValueError("BBR025: excitation.normal_entry_theta_deg must be 180")
    _expect(_req(exc, "polarization_convention", "excitation.", "str"), POLARIZATION,
            "excitation.polarization_convention")
    for key in ("incident_phi_deg", "incident_theta_deg"):
        v = _req(exc, key, "excitation.", "list")
        if not v or not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
            raise ValueError(f"BBR024: excitation.{key} must be a non-empty list of numbers")

    ff = _req(sc, "far_field", "", "dict")                                         # F6
    _expect(_req(ff, "coordinate_system", "far_field.", "str"), "exit_cs", "far_field.coordinate_system")
    _expect(_req(ff, "definition", "far_field.", "str"), "Theta-Phi", "far_field.definition")
    _expect(_req(ff, "component_basis", "far_field.", "str"), "spherical_in_exit_cs",
            "far_field.component_basis")
    _axis_range(ff, "theta_deg", "far_field.")
    _axis_range(ff, "phi_deg", "far_field.")
    _req(ff, "points_per_key", "far_field.", "int")
    if _req(ff, "columns", "far_field.", "list") != FF_COLUMNS:
        raise ValueError("BBR025: far_field.columns differ from the positional order BBRHFSSData reads")

    xf = _req(sc, "exit_field", "", "dict")                                        # F7
    _expect(_req(xf, "coordinate_system", "exit_field.", "str"), "exit_cs", "exit_field.coordinate_system")
    if not _req(xf, "points_in_si", "exit_field.", "bool"):
        raise ValueError("BBR025: exit_field.points_in_si must be true")
    _expect(_req(xf, "plane", "exit_field.", "str"), "x_e=0", "exit_field.plane")
    in_ref = _req(xf, "field_in_ref_cs", "exit_field.", "bool")
    if _req(xf, "field_components_frame", "exit_field.", "str") != ("exit_cs" if in_ref else "hfss_global"):
        raise ValueError("BBR025: exit_field.field_in_ref_cs and field_components_frame contradict each other")
    grid = _req(xf, "grid", "exit_field.", "dict")
    _axis_range(grid, "y_e", "exit_field.grid.")
    _axis_range(grid, "z_e", "exit_field.grid.")
    _req(xf, "points_per_key_retained", "exit_field.", "int")
    if _req(xf, "outside_points", "exit_field.", "str") not in ("none", "omitted", "zero"):
        raise ValueError("BBR024: exit_field.outside_points must be none, omitted or zero")
    if _req(xf, "columns", "exit_field.", "list") != (WG_COLUMNS_2 if major == 2 else WG_COLUMNS):
        raise ValueError("BBR025: exit_field.columns differ from the positional order BBRHFSSData reads")

    cs = _req(xf, "cross_section", "exit_field.", "dict")                         # F8
    shape = _req(cs, "shape", "exit_field.cross_section.", "str")
    if shape == "rectangle":
        for k in ("y_e_half_m", "z_e_half_m"):
            if not _req(cs, k, "exit_field.cross_section.", "num") > 0:
                raise ValueError(f"BBR024: exit_field.cross_section.{k} must be > 0")
    elif shape == "disc":
        if not _req(cs, "radius_m", "exit_field.cross_section.", "num") > 0:
            raise ValueError("BBR024: exit_field.cross_section.radius_m must be > 0")
    elif shape == "polygon":
        raise ValueError("BBR025: cross_section shape polygon is reserved; BBRsim does not support it yet")
    else:
        raise ValueError(f"BBR024: unknown cross_section shape {shape!r}")

    if major == 1:                                                                 # F9
        tr = _req(sc, "transmittance", "", "dict")
        _expect(_req(tr, "definition", "transmittance.", "str"), TRANSMITTANCE, "transmittance.definition")
        if _req(tr, "incoming_includes_cos_theta", "transmittance.", "bool"):
            raise ValueError("BBR025: transmittance.incoming_includes_cos_theta must be false")
    else:                      # under 2.x the transmittance, exit-field, cross-term and per-key blocks replace it
        _schema2(sc, cs)

    sym = _req(sc, "symmetry", "", "dict")                                         # F10
    for k in ("mirror_l", "mirror_g", "end_to_end"):
        if not _req(sym, k, "symmetry.", "bool"):
            raise ValueError(f"BBR025: symmetry.{k} is false; the sampler folds by both mirrors "
                             "and serves both ends from one table")
    _req(sym, "rotational", "symmetry.", "bool")
    # the boundary model the transport implements: radiation end faces, PEC walls; non-transmission is
    # specular reflection with no wall loss, and the mode formulas assume a closed PEC guide (CDX-003)
    bd = _req(sc, "boundaries", "", "dict")
    _expect(_req(bd, "entrance", "boundaries.", "str"), "radiation", "boundaries.entrance")
    _expect(_req(bd, "exit", "boundaries.", "str"), "radiation", "boundaries.exit")
    _expect(_req(bd, "walls", "boundaries.", "str"), "PEC", "boundaries.walls")

    # the gap must be vacuum (CDX-029); null or an absent key is accepted (the legacy sidecars carry null)
    prov = sc.get("provenance")
    if isinstance(prov, dict) and "material" in prov:
        mat = prov["material"]
        if isinstance(mat, str):
            if mat.lower() != "vacuum":
                raise ValueError(f"BBR025: provenance.material is {mat!r}; the gap must be vacuum")
        elif mat is not None:
            raise ValueError(f"BBR024: provenance.material must be a string or null, not {type(mat).__name__}")

    ext = _req(_req(sc, "geometry", "", "dict"), "extent_mm", "geometry.", "dict")  # F11
    for k in "plg":
        if not _req(ext, k, "geometry.extent_mm.", "num") > 0:
            raise ValueError(f"BBR024: geometry.extent_mm.{k} must be > 0")
    # a rectangle's l and g are the full widths of its declared section (CDX-028); a disc's are not
    # compared with its diameter (the round gap's extent_mm is the HFSS bounding box); p is F12's, C++ only
    if cs["shape"] == "rectangle":
        for k, half in (("l", "y_e_half_m"), ("g", "z_e_half_m")):
            width = 2e3 * cs[half]
            if abs(ext[k] - width) > 1e-6 * width:
                raise ValueError(f"BBR024: geometry.extent_mm.{k} = {ext[k]!r} mm differs from the full width "
                                 f"{width!r} mm of exit_field.cross_section (2 x {half}) by more than 1e-6 relative")

    md = _req(sc, "modes", "", "dict")                                             # F13
    got = (_req(md, "mode", "modes.", "str"), _req(md, "cutoff_ghz", "modes.", "num"),
           _req(md, "polarization_filter_limit_ghz", "modes.", "num"))
    want = lowest_mode(cs)
    if got[0] != want[0] or abs(got[1] - want[1]) > MODE_REL * want[1] or abs(got[2] - want[2]) > MODE_REL * want[2]:
        raise ValueError(f"BBR025: modes (mode, cutoff_ghz, polarization_filter_limit_ghz) = {got}, "
                         f"the declared cross-section gives {want}")
    if major == 2:             # OR-031: the declared cutoffs, re-derived entry by entry
        _check_regime_cutoffs(sc, md, cs, got[2])


def _contains(section, y, z):
    if section["shape"] == "rectangle":
        return (abs(y) <= section["y_e_half_m"] * (1 + SECTION_REL)
                and abs(z) <= section["z_e_half_m"] * (1 + SECTION_REL))
    return y * y + z * z <= section["radius_m"] ** 2 * (1 + SECTION_REL)


def _nodes(values, axis, tol):
    """The node index of each value on a declared axis, node i at min + i (max - min)/(count - 1),
    or -1 off it: matched to 1e-6 of the step, to tol when the axis has one node (as BBRHFSSData)."""
    v = np.asarray(values, float)
    lo, hi, n = axis["min"], axis["max"], axis["count"]
    if n <= 1:
        return np.where(np.abs(v - lo) <= tol, 0, -1)
    step = (hi - lo) / (n - 1)
    i = np.rint((v - lo) / step).astype(int)
    return np.where((i >= 0) & (i < n) & (np.abs(v - (lo + i * step)) <= 1e-6 * step), i, -1)


def _node_at(axis, i):
    return axis["min"] if axis["count"] <= 1 else axis["min"] + i * (axis["max"] - axis["min"]) / (axis["count"] - 1)


def _axis_text(name, axis):
    return f"{name} {axis['min']!r} to {axis['max']!r} in {axis['count']} nodes"


def _check_freq(stem, d):
    """Every Freq token parses and agrees with the directory frequency to 0.1 %, as BBRHFSSData checks
    each row (BBR009; CDX-001)."""
    tag = f"{stem}_Ephi={d['ephi']}"
    f_dir = hfss.parse_frequency_GHz(stem.rsplit("_", 1)[1])
    for name, tokens in d["freq_tokens"].items():
        for tok in tokens:
            f = hfss.parse_frequency_GHz(tok)
            if f is None or f_dir is None or abs(f - f_dir) > 1e-3 * f_dir:
                raise ValueError(f"BBR009: {tag}/{name}: Freq column {tok!r} does not match the directory "
                                 f"frequency {f_dir} GHz (tolerance 0.1 %)")


def _far_rel(a, b):
    """a and b differ by more than 1e-6 of the larger magnitude (BBRHFSSData's FarRel)."""
    return abs(a - b) > DIVISOR_REL * max(abs(a), abs(b))


def _check_rows_v2(stem, d, sc, recs, sidecar_path):
    """The row checks of BBRHFSSData::LoadWaveguide under 2.x (BBR012, naming key, polarization and value),
    in file order: InterceptFactor constant within its key, equal to the key's c_key and to the closed form
    of the declared grid; OutgoingPower exactly 0 on every row of a dark basis (OR-038); and T = OutgoingPower
    / (IngoingPower * InterceptFactor) from the key's first row against the sidecar's t_flux / c_k, 1e-6
    relative on a lit key, both exactly 0 on a dark one. A key without a record is left to C2."""
    ephi = d["ephi"]
    thetas = [float(x) for x in sc["excitation"]["incident_theta_deg"]]
    iw, out, inp, fac = d["rows2"]

    def find(phi, theta):
        for r in recs:
            if abs(r["phi"] - phi) <= KEY_TOL and abs(r["theta"] - theta) <= KEY_TOL:
                return r
        return None

    first, cache = {}, {}
    for n in range(len(out)):
        phi, theta = float(iw[n, 0]), float(iw[n, 1])
        key = (round(phi, 2), round(theta, 2))
        if (phi, theta) not in cache:
            cache[(phi, theta)] = find(phi, theta)
        rec = cache[(phi, theta)]
        f, o, i = float(fac[n]), float(out[n]), float(inp[n])

        def bad(what):
            raise ValueError(f"BBR012: {stem}_Ephi={ephi}/waveguide.csv, key ({_rep(phi)}, {_rep(theta)}) "
                             f"polarization {ephi}: {what} (sidecar {sidecar_path})")

        if key in first and _far_rel(f, first[key]):
            bad(f"InterceptFactor {_rep(f)} differs from the key's first row, {_rep(first[key])}, by more than "
                "1e-6 relative: not constant within the key")
        dark = False
        if rec is not None:
            if _far_rel(f, rec["c_k"]):
                bad(f"InterceptFactor {_rep(f)} differs from the sidecar's transmittance.c_key {_rep(rec['c_k'])} "
                    "by more than 1e-6 relative")
            closed = intercept_factor(rec["theta"], thetas)
            if _far_rel(f, closed):
                bad(f"InterceptFactor {_rep(f)} differs from {_rep(closed)}, the closed form of the declared grid, "
                    "by more than 1e-6 relative")
            dark = rec["dark1"] if ephi else rec["dark0"]
            if dark and o != 0.0:
                bad(f"the sidecar declares this basis dark, but OutgoingPower = {_rep(o)}; a dark basis has "
                    "OutgoingPower exactly 0 on every row of its key (OR-038)")
        if key not in first:
            t = o / (i * f) if i * f > 0.0 else 0.0
            first[key] = f
            if rec is not None:
                tf = (rec["t_flux1"] if ephi else rec["t_flux0"]) / rec["c_k"]
                if (t != 0.0 or tf != 0.0) if dark else _far_rel(t, tf):
                    bad(f"T = OutgoingPower / (IngoingPower * InterceptFactor) = {_rep(t)} but the sidecar's "
                        f"t_flux / c_k = {_rep(tf)}" + ("; a dark basis has both exactly 0" if dark
                                                         else ", more than 1e-6 relative apart"))


def _check_one(stem, d, sc, exempt):
    """C1-C5 on the description of one polarization; exempt holds the (key, Y, Z) that C5 lets
    lie outside the section."""
    tag = f"{stem}_Ephi={d['ephi']}"
    ff, xf = sc["far_field"], sc["exit_field"]
    if d["far_field"]["columns"] != ff["columns"] or d["exit_field"]["columns"] != xf["columns"]:
        raise ValueError(f"BBR013: {tag}: a CSV header differs from the sidecar's columns")              # C1
    _check_freq(stem, d)
    # C2: keys rounded to 0.01 deg, as BBRHFSSData rounds them; the waveguide keys must equal the
    # declared grid, and the far-field keys must equal them (BBR007 in BBRHFSSData too). A declared
    # key without rows would be served by its nearest neighbour, another incidence's tables (CDX-021).
    phis = {round(float(v), 2) for v in sc["excitation"]["incident_phi_deg"]}
    thetas = {round(float(v), 2) for v in sc["excitation"]["incident_theta_deg"]}
    off = sorted(k for k in d["keys"] if k[0] not in phis or k[1] not in thetas)
    if off:
        raise ValueError(f"BBR007: {tag}: keys {off} are not on the declared incidence grid")
    if d["ff_keys"] != d["keys"]:
        raise ValueError(f"BBR007: {tag}: far_field.csv and waveguide.csv disagree on incidence keys "
                         f"(far field only {sorted(d['ff_keys'] - d['keys'])}, "
                         f"waveguide only {sorted(d['keys'] - d['ff_keys'])})")
    missing = sorted({(p, t) for p in phis for t in thetas} - d["keys"])
    if missing:
        raise ValueError(f"BBR007: {tag}: declared incidence keys {missing} have no rows; the nearest-key "
                         f"lookup would serve another incidence's tables for them")
    for name in ("theta_deg", "phi_deg"):                                                               # C3
        got, want = d["far_field"][name], ff[name]
        if got["count"] != want["count"] or got["min"] < want["min"] - 1e-9 or got["max"] > want["max"] + 1e-9:
            raise ValueError(f"BBR012: {tag}: far-field {name} {got} outside or unlike the declared {want}")
        if got["count"] > 1:                     # the step, which only this validator checks
            step = _req(want, "step", f"far_field.{name}.", "num")
            if abs(step - got["step"]) > 1e-9 * got["step"]:
                raise ValueError(f"BBR012: {tag}: far-field {name} step {step} declared, the data give "
                                 f"{got['step']} = (max - min)/(count - 1)")
    if d["far_field"]["points_per_key"] != ff["points_per_key"]:
        raise ValueError(f"BBR012: {tag}: {d['far_field']['points_per_key']} far-field rows per key, "
                         f"declared {ff['points_per_key']}")
    # C4: with "none" and "zero" every lattice point is written, so the retained distinct counts
    # equal the declared ones; with "omitted" the producer declares the whole export lattice and
    # drops the points outside the section, so a round face whose radius is not a multiple of the
    # step leaves its outermost lattice columns without a retained point: the retained counts may
    # then be smaller than the declared ones, never larger (BBRHFSSData::CheckAgainstSidecar).
    omitted = xf["outside_points"] == "omitted"
    for axis in ("y_e", "z_e"):                                                                         # C4
        got, want = d["exit_field"]["grid"][axis], xf["grid"][axis]
        tol = 1e-9 * max(abs(want["min"]), abs(want["max"]))
        if (got["min"] < want["min"] - tol or got["max"] > want["max"] + tol
                or (not omitted and got["count"] != want["count"])):
            raise ValueError(f"BBR012: {tag}: exit grid {axis} {got} outside or unlike the declared {want}")
        if omitted and got["count"] > want["count"]:
            raise ValueError(f"BBR012: {tag}: exit grid {axis} {got} has more distinct values than the declared "
                             f"lattice {want} (with outside_points omitted the retained count may be smaller, "
                             f"not larger)")
    if d["exit_field"]["points_per_key_retained"] != xf["points_per_key_retained"]:
        raise ValueError(f"BBR012: {tag}: exit points per key differ from points_per_key_retained")
    if d["x_max_abs_m"] > 1e-12:
        raise ValueError(f"BBR012: {tag}: exit points off the plane x_e = 0 (|X| up to {d['x_max_abs_m']})")
    outside = [(y, z) for k, (y, z) in zip(d["rows_key"], d["yz_m"])                                   # C5
               if not _contains(xf["cross_section"], y, z) and (k, float(y), float(z)) not in exempt]
    if outside:
        raise ValueError(f"BBR025: {tag}: {len(outside)} exit point(s) outside the declared "
                         f"cross-section, e.g. {outside[0]}")
    # C3, C4 per key (CDX-022): every far-field row on a node of the declared Phi x Theta grid and
    # every exit point on a node of the declared Y x Z lattice, no node twice in one key, every key
    # holding all of them (the whole exit lattice, or with "omitted" its nodes inside the section).
    # The dataset-wide checks above cannot see a repeated row, an uneven axis or keys that differ
    # in compensating ways.
    phi_ax, theta_ax = ff["phi_deg"], ff["theta_deg"]
    y_ax, z_ax = xf["grid"]["y_e"], xf["grid"]["z_e"]
    y_tol = 1e-9 * max(abs(y_ax["min"]), abs(y_ax["max"]))
    z_tol = 1e-9 * max(abs(z_ax["min"]), abs(z_ax["max"]))
    lattice = {(i, j) for i in range(y_ax["count"]) for j in range(z_ax["count"])
               if not omitted or _contains(xf["cross_section"], _node_at(y_ax, i), _node_at(z_ax, j))}
    for k, rows in d["ff_rows"].items():
        i, j = _nodes(rows[:, 0], phi_ax, 1e-9), _nodes(rows[:, 1], theta_ax, 1e-9)
        off = np.nonzero((i < 0) | (j < 0))[0]
        if len(off):
            r = rows[off[0]]
            raise ValueError(f"BBR012: {tag}: far-field row (Phi, Theta) = ({float(r[0])!r}, {float(r[1])!r}) of key {k} is not a "
                             f"node of the declared grid {_axis_text('Phi', phi_ax)} x {_axis_text('Theta', theta_ax)}")
        seen = set()
        for n, pair in enumerate(zip(i.tolist(), j.tolist())):
            if pair in seen:
                raise ValueError(f"BBR012: {tag}: far-field point (Phi, Theta) = ({float(rows[n, 0])!r}, {float(rows[n, 1])!r}) "
                                 f"appears twice in key {k}")
            seen.add(pair)
        if len(seen) != phi_ax["count"] * theta_ax["count"]:
            raise ValueError(f"BBR012: {tag}: key {k} holds {len(seen)} of the {phi_ax['count'] * theta_ax['count']} "
                             f"nodes of the declared far-field grid")
    for k, rows in d["exit_rows"].items():
        i, j = _nodes(rows[:, 1], y_ax, y_tol), _nodes(rows[:, 2], z_ax, z_tol)
        off = np.nonzero((i < 0) | (j < 0))[0]
        if len(off):
            r = rows[off[0]]
            raise ValueError(f"BBR012: {tag}: exit point (Y, Z) = ({float(r[1])!r}, {float(r[2])!r}) m of key {k} is not a node of "
                             f"the declared lattice {_axis_text('Y', y_ax)} x {_axis_text('Z', z_ax)}")
        seen = set()
        for n, pair in enumerate(zip(i.tolist(), j.tolist())):
            if pair in seen:
                raise ValueError(f"BBR012: {tag}: exit point (Y, Z) = ({float(rows[n, 1])!r}, {float(rows[n, 2])!r}) m appears "
                                 f"twice in key {k}")
            seen.add(pair)
        if seen != lattice:
            raise ValueError(f"BBR012: {tag}: key {k} holds {len(seen)} exit points; every key holds each of the "
                             f"{len(lattice)} declared lattice nodes{' inside the section' if omitted else ''} once")


def check_csvs(base, stem, sc):
    """C1-C5 on both polarizations, as BBRHFSSData runs them when it is given the sidecar, plus
    the far-field step and the agreement of Ephi=1 with Ephi=0 on keys, grids and per-key counts. Under
    2.x the tables have 16 columns and the row checks of the loader come first (_check_rows_v2).

    Precondition: sc has passed check(); its blocks and types are not re-checked here. A missing
    far_field.csv or waveguide.csv raises ValueError starting BBR001 or BBR002, as in BBRHFSSData.
    """
    v2 = schema_major(sc) == 2
    d0, d1 = describe_csvs(base, stem, v2)["polarizations"]
    if v2:               # the row checks of the 2.x loader come before the declaration checks, per polarization
        recs = key_records(sc)
        for p in (d0, d1):
            _check_freq(stem, p)
            _check_rows_v2(stem, p, sc, recs, path_for(base, stem))
    # C5: with outside_points "zero" a point outside the section is allowed when its six field
    # components are exactly zero in both polarizations; with "none" or "omitted" none is.
    exempt = d0["zero_points"] & d1["zero_points"] if sc["exit_field"]["outside_points"] == "zero" else set()
    for p in (d0, d1):
        _check_one(stem, p, sc, exempt)
    agree = (("incidence keys", lambda p: p["keys"]),
             ("far-field Theta", lambda p: p["far_field"]["theta_deg"]),
             ("far-field Phi", lambda p: p["far_field"]["phi_deg"]),
             ("far-field rows per key", lambda p: p["far_field"]["points_per_key"]),
             ("exit grid", lambda p: p["exit_field"]["grid"]),
             ("exit points per key", lambda p: p["exit_field"]["points_per_key_retained"]))
    for what, get in agree:
        if get(d0) != get(d1):
            raise ValueError(f"BBR012: {stem}: Ephi=1 {what} {get(d1)} differ from Ephi=0's {get(d0)}")
    # BBRHFSSData pairs each key's Ephi=1 rows with its Ephi=0 rows by order and requires each pair at
    # the same far-field (Phi, Theta) to 1e-9 deg and exit (X, Y, Z) to 1e-12 m (BBR012). A polarization
    # with its rows reordered passes every check above (CDX-002).
    for what, rows, names, tol in (("far-field row", "ff_rows", "(Phi, Theta)", 1e-9),
                                   ("exit point", "exit_rows", "(X, Y, Z)", 1e-12)):
        for k, a in d0[rows].items():
            b = d1[rows][k]
            bad = np.nonzero(np.any(np.abs(a - b) > tol, axis=1))[0]
            if len(bad):
                n = int(bad[0])
                raise ValueError(f"BBR012: {stem}: key {k}: Ephi=1 {what} {n} lies at {names} = "
                                 f"{tuple(float(v) for v in b[n])}, the Ephi=0 one at {tuple(float(v) for v in a[n])}; "
                                 f"BBRHFSSData pairs the two polarizations' rows by order")


def _same(got, want, where="modes"):
    """got carries every key of want with an equal value (floats to MODE_REL); extra keys, which a
    later 1.x schema may add, are ignored."""
    if isinstance(want, dict):
        if not isinstance(got, dict) or not set(want) <= set(got):
            raise ValueError(f"BBR025: {where}: keys {sorted(got) if isinstance(got, dict) else got} "
                             f"lack some of {sorted(want)}")
        for k in want:
            _same(got[k], want[k], f"{where}.{k}")
    elif isinstance(want, list):
        if not isinstance(got, list) or len(got) != len(want):
            raise ValueError(f"BBR025: {where}: {len(got) if isinstance(got, list) else got} entries, "
                             f"expected {len(want)}")
        for i, (g, w) in enumerate(zip(got, want)):
            _same(g, w, f"{where}[{i}]")
    elif isinstance(want, float):
        if not isinstance(got, (int, float)) or abs(got - want) > MODE_REL * max(abs(want), 1e-300):
            raise ValueError(f"BBR025: {where} = {got}, expected {want}")
    elif got != want:
        raise ValueError(f"BBR025: {where} = {got!r}, expected {want!r}")


def check_full(base, dataset_id, stem, dir_freq_ghz):
    """Everything: check(), under 1.x the full modes block (basis is recorded only, not compared), the
    sha256, bytes and rows of the four CSVs, and check_csvs(). Unknown fields are ignored. Under 2.x the
    modes block is the regime-cutoff list check() re-derived."""
    sc = load(base, stem)
    check(sc, dataset_id, stem, dir_freq_ghz)
    major = schema_major(sc)
    if major == 1:
        want = modes_for(sc["exit_field"]["cross_section"], sc["frequency_ghz"])
        _same(sc["modes"], {k: v for k, v in want.items() if k != "basis"})
    files = _req(sc, "files", "", "dict")
    for rel, entry in file_entries(base, stem).items():
        got = files.get(rel)
        if not isinstance(got, dict):
            raise ValueError(f"BBR024: {path_for(base, stem)}: files has no entry for {rel}")
        for k in ("sha256", "bytes", "rows"):
            if got.get(k) != entry[k]:
                raise ValueError(f"BBR024: {path_for(base, stem)}: files[{rel!r}].{k} is {got.get(k)!r}, "
                                 f"the CSV gives {entry[k]!r}")
    check_csvs(base, stem, sc)
    # the load-time checks of the BBRHFSSData constructor, BBR027 included (CDX-008): a dataset that passes
    # here but refuses to load would stop a Geant4 run at the first photon of its frequency. The Python
    # mirror reads both schemas: under 2.x by the 2.0 positions, with the weights, the per-photon cap and
    # the dark masking (OR-034, OR-038, HF-040, HF-042).
    hfss.load_dataset(stem, base)
    return sc


# --- the schema-2 builder ------------------------------------------------------------------------

def _hemisphere_power_w(rows, dtheta_deg, dphi_deg):
    """Far-field power over the outward hemisphere the CSV spans: rows with sin(Theta) cos(Phi) > 1e-6
    summed with weight sin(Theta) dTheta dPhi (radians) and (|rE_theta|^2 + |rE_phi|^2) / (2 eta0)."""
    theta, phi = np.radians(rows.Theta.to_numpy(float)), np.radians(rows.Phi.to_numpy(float))
    sin_t = np.sin(theta)
    outward = sin_t * np.cos(phi) > 1e-6
    inten = (rows.rEtheta_real.to_numpy(float) ** 2 + rows.rEtheta_imag.to_numpy(float) ** 2
             + rows.rEphi_real.to_numpy(float) ** 2 + rows.rEphi_imag.to_numpy(float) ** 2)
    weight = sin_t * math.radians(dtheta_deg) * math.radians(dphi_deg)
    return float(np.sum((weight * inten)[outward]) / (2.0 * ETA0))


def build_schema2_from_csvs(base, dataset_id, frequency_label, *, section, extent_mm, provenance,
                            exit_origin_mm, plane_wave_origin_mm, bounding_box_mm, step_m=None, w=None,
                            rotational=False, abs_powers=None, pose_rule="canonical-z", resolution_mm=None,
                            outside_points="none", rim_points="included", bounds_method="vertices", walls="PEC"):
    """A schema-2.0 sidecar for <dataset_id>_<frequency_label>, read off its finalized 2.0 CSVs.

    c_key and cell_90 come from the declared grid (the incident angles the tables hold); t_flux per key
    and polarization is the first row's OutgoingPower over its IngoingPower (the magnitude numerator is
    copied from OutgoingPower unless abs_powers {(phi, theta): (p0, p1)} gives it); W = 1 unless w
    {(phi, theta): (w0, w1)} gives it, and always 1 on a dark basis; ff_over_flux is the outward-hemisphere
    far-field power over the stored exit flux, null on a dark basis. step_m: a number (both axes), a dict
    {"y_e": ..., "z_e": ...}, or None to read it off the lattice, (max - min)/(count - 1).

    A basis is dark where |t_flux| <= 1e-12 (the producer's threshold). Its OutgoingPower must be exactly 0
    on every row: numerical residue is refused (BBR012) rather than written into a sidecar that contradicts
    its own tables; the caller converts a numerical zero before building. The CSVs are never modified, and
    the checksums are those of the files as they stand.
    """
    stem = f"{dataset_id}_{frequency_label}"
    freq = hfss.parse_frequency_GHz(frequency_label)
    d = describe_csvs(base, stem, True)
    d0, d1 = d["polarizations"]
    phis, thetas = d["incident_phi_deg"], d["incident_theta_deg"]
    cell = cell_90(thetas)
    if cell is None:
        raise ValueError("BBR025: the declared incident_theta_deg has no theta = 90 deg entry with a theta above "
                         "it: cell_90 cannot be derived (HF-049)")
    delta = _delta(thetas)
    mode, cut, filt = lowest_mode(section)
    far = [_read_csv(base, stem, e, "far_field.csv", FF_COLUMNS) for e in (0, 1)]
    far_groups = [f.assign(kp=f.IWavePhi.round(2), kt=f.IWaveTheta.round(2)).groupby(["kp", "kt"]) for f in far]
    rows = []
    for dd in (d0, d1):
        iw = dd["rows2"][0]
        keys = pd.DataFrame({"kp": np.round(iw[:, 0], 2), "kt": np.round(iw[:, 1], 2)})
        rows.append(keys.groupby(["kp", "kt"]).indices)
    c_key, per_key, scale, ff_ratio = [], [], [], []
    for ph in phis:
        for t in thetas:
            key = (round(ph, 2), round(t, 2))
            ck = intercept_factor(t, thetas)
            tflux, absp, darks, ratios = [], [], [], []
            for e, dd in enumerate((d0, d1)):
                idx = rows[e].get(key)
                if idx is None:
                    raise ValueError(f"BBR007: {stem}_Ephi={e}: declared key ({_rep(ph)}, {_rep(t)}) has no rows")
                out, inp = dd["rows2"][1], dd["rows2"][2]
                flux, a_i = float(out[idx[0]]), float(inp[idx[0]])
                if not a_i > 0.0:
                    raise ValueError(f"BBR012: {stem}_Ephi={e}: key ({_rep(ph)}, {_rep(t)}) has IngoingPower "
                                     f"{_rep(a_i)}; the entrance power must be positive")
                tf = flux / a_i
                dark = abs(tf) <= DARK_REL
                if dark and np.any(out[idx] != 0.0):
                    raise ValueError(f"BBR012: {stem}_Ephi={e}: key ({_rep(ph)}, {_rep(t)}) polarization {e} is dark "
                                     f"(|t_flux| = {_rep(abs(tf))} <= 1e-12) but holds OutgoingPower other than "
                                     "exactly 0 (e.g. " + _rep(float(out[idx][out[idx] != 0.0][0])) + "); convert "
                                     "the numerical zero before building (OR-038)")
                tflux.append(0.0 if dark else tf)
                absp.append(float(abs_powers[(ph, t)][e]) if abs_powers and (ph, t) in abs_powers else flux)
                darks.append(dark)
                ratios.append(None if dark else _hemisphere_power_w(
                    far_groups[e].get_group(key), d["far_field"]["theta_deg"]["step"], d["far_field"]["phi_deg"]["step"]) / flux)
            ww = w.get((ph, t)) if w else None
            w0, w1 = (1.0 if darks[0] else (ww[0] if ww else 1.0)), (1.0 if darks[1] else (ww[1] if ww else 1.0))
            c_key.append({"phi_deg": ph, "theta_deg": t, "c_k": ck})
            per_key.append({"phi_deg": ph, "theta_deg": t, "c_k": ck, "polarizations": {
                "0": {"t_flux": tflux[0], "outgoing_power_abs_re_s_w": absp[0]},
                "1": {"t_flux": tflux[1], "outgoing_power_abs_re_s_w": absp[1]}}})
            scale.append({"phi_deg": ph, "theta_deg": t, "w0": w0, "w1": w1, "dark0": darks[0], "dark1": darks[1]})
            ff_ratio.append({"phi_deg": ph, "theta_deg": t, "polarizations": {"0": ratios[0], "1": ratios[1]}})
    grid = {}
    for ax in ("y_e", "z_e"):
        g = dict(d["exit_field"]["grid"][ax])
        step = step_m.get(ax) if isinstance(step_m, dict) else step_m
        if step is None:
            if g["count"] < 2:
                raise ValueError(f"BBR024: exit_field.grid.{ax} has {g['count']} node; its step cannot be derived")
            step = (g["max"] - g["min"]) / (g["count"] - 1)
        g["step"] = float(step)
        grid[ax] = g
    return {
        "schema_version": "2.0",
        "dataset_id": dataset_id,
        "frequency_ghz": freq,
        "frequency_label": frequency_label,
        "provenance": provenance,
        "model": {"kind": SCHEMA2_STRINGS["model.kind"]},
        "conventions": {"time_dependence": "exp(+j omega t)", "phasor": "peak"},
        "frames": {
            "canonical": CANONICAL, "pose_rule": pose_rule,
            "hfss_global_axes_in_canonical": HFSS_AXES,
            "entrance_face": {"axis": "z", "side": "min"}, "exit_face": {"axis": "z", "side": "max"},
            "entrance_outward_normal_global": [0, 0, -1], "exit_outward_normal_global": [0, 0, 1],
            "exit_cs": {"name": "outgoing_cs", "origin": "exit_face_center",
                        "origin_mm_global": exit_origin_mm, **EXIT_CS_GLOBAL},
            "exit_cs_axes_in_canonical": IDENTITY,
        },
        "excitation": {
            "coordinate_system": "global", "incidence_convention": "arrival_direction",
            "normal_entry_theta_deg": 180.0,
            "incident_phi_deg": phis, "incident_theta_deg": thetas,
            "ei_v_per_m": 1.0, "incoming_power_w": d["incoming_power_w"],
            "polarization_convention": POLARIZATION,
            "plane_wave_origin": "entrance_face_center", "origin_mm_global": plane_wave_origin_mm,
        },
        "far_field": {"coordinate_system": "exit_cs", "definition": "Theta-Phi",
                      "component_basis": "spherical_in_exit_cs", "radiation_surface": "exit_face",
                      **d["far_field"], "normalization": SCHEMA2_STRINGS["far_field.normalization"],
                      "ff_over_flux": ff_ratio},
        "exit_field": {
            "coordinate_system": "exit_cs", "points_in_si": True, "field_in_ref_cs": False,
            "field_components_frame": "hfss_global", "plane": "x_e=0", "resolution_mm": resolution_mm,
            "grid": grid, "cross_section": section, "outside_points": outside_points,
            "bounds_method": bounds_method, "rim_points": rim_points,
            "points_per_key_retained": d["exit_field"]["points_per_key_retained"],
            "columns": d["exit_field"]["columns"],
            "tangential_components": list(SCHEMA2_TANGENTIAL),
            "quadrature": SCHEMA2_QUADRATURE[section["shape"]],
            "sampler": SCHEMA2_STRINGS["exit_field.sampler"],
            "scale": {"rule": SCHEMA2_STRINGS["exit_field.scale.rule"], "per_key": scale},
            "column_order_check": {"rule": SCHEMA2_STRINGS["exit_field.column_order_check.rule"],
                                   "method": "reference",
                                   "reference": {"dataset_id": "InfParallelPlate_crack1Rohan",
                                                 "frequency_ghz": 500, "passed": True}},
            "plane_tolerance_m": SCHEMA2_PLANE_TOL_M,
        },
        "transmittance": {
            "definition": SCHEMA2_STRINGS["transmittance.definition"],
            "numerator": SCHEMA2_STRINGS["transmittance.numerator"],
            "incoming": SCHEMA2_STRINGS["transmittance.incoming"],
            "above_one": SCHEMA2_STRINGS["transmittance.above_one"],
            "c_key": c_key, "cell_90": cell,
            "cell_90_assumes": {"emitter": "uniform_theta", "cos_floor": COS_FLOOR,
                                "band_deg": [90.0 - delta / 2, 90.0]},
            "per_key": per_key,
        },
        "cross_term": {"density": SCHEMA2_STRINGS["cross_term.density"], "definition": SCHEMA2_CROSS_TERM,
                       "valid_below_ghz": filt, "above_limit": SCHEMA2_STRINGS["cross_term.above_limit"]},
        "symmetry": {"mirror_l": True, "mirror_g": True, "end_to_end": True, "rotational": rotational},
        "boundaries": {"entrance": "radiation", "exit": "radiation", "walls": walls},
        "geometry": {"shape": "box" if section["shape"] == "rectangle" else "cylinder",
                     "extent_mm": extent_mm, "bounding_box_mm": bounding_box_mm},
        "modes": {"mode": mode, "cutoff_ghz": cut, "polarization_filter_limit_ghz": filt,
                  "regime_cutoffs_ghz": regime_cutoffs(section)},
        "files": file_entries(base, stem),
    }
