"""bbrsim.sidecar - the HFSS dataset sidecar <stem>.dataset.json (schema 1.x).

Python twin of library/src/BBRDatasetSidecar.cc (checks F1-F11 and F13, codes
BBR024 and BBR025) and of the sidecar checks in BBRHFSSData (C1-C5, on both
polarizations), plus what the C++ does not do: the full mode lists (Bessel zeros
through scipy), the CSV checksums, the far-field step, the agreement of the two
polarizations, and the builders the mock generators and the legacy-sidecar
script use. The schema (agreed with Blackbody-Simulations on 2026-10-05) is described
in validation/README.md, Dataset sidecars. Errors are ValueError whose message
starts with the C++ code. invariant() and same_invariant() hold what every
frequency of one ID must share: physics only, numbers to 1e-9 relative, over the
fields both sidecars carry; invariant_diff() names the blocks that differ.

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


def disc_modes(radius_m, freq_ghz, limit_ghz=LIST_LIMIT_GHZ):
    """Closed ideal-PEC circular guide: TE_np = x'_np c/(2 pi R), TM_np = x_np c/(2 pi R)."""
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
        if not added:          # the first zero grows with m, so no higher m qualifies
            break
        m += 1
    cutoffs = [{"mode": mode_label(k, m, p), "cutoff_ghz": f, "degeneracy": 1 if m == 0 else 2}
               for f, k, m, p in _order(entries)]
    lowest = X11_PRIME * scale
    return {"cutoff_ghz": lowest, "mode": "TE11", "list_limit_ghz": limit_ghz,
            "polarization_filter_limit_ghz": lowest, "cutoffs": cutoffs,
            "propagating_count": sum(1 for c in cutoffs if c["cutoff_ghz"] <= freq_ghz)}


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
    if df.isna().to_numpy().any():
        raise ValueError(f"BBR013: {p}: a row lacks fields or holds an empty field")
    return df


def _describe_one(base, stem, ephi):
    """What the two CSVs of one polarization contain: header, incidence keys, grids, per-key counts."""
    ff = _read_csv(base, stem, ephi, "far_field.csv", FF_COLUMNS)
    wg = _read_csv(base, stem, ephi, "waveguide.csv", WG_COLUMNS)

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
                        zip(rows_key, yz, (wg[FIELD_COLUMNS].to_numpy(float) == 0).all(axis=1)) if zero},
        "rows_key": rows_key,
        "keys": set(rows_key),
        "ff_keys": {(round(float(p), 2), round(float(t), 2)) for p, t in
                    ff[["IWavePhi", "IWaveTheta"]].itertuples(index=False)},
    }


def describe_csvs(base, stem):
    """What the CSVs of <stem> contain. The top-level fields describe Ephi=0 (what the builder
    records); "polarizations" holds the description of each of Ephi=0 and Ephi=1."""
    d0, d1 = _describe_one(base, stem, 0), _describe_one(base, stem, 1)
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
    each file: the frame mapping, symmetry, boundaries, geometry shape and extents,
    the modes block without basis (recorded only) and propagating_count (per
    frequency), and the exit cross-section. Descriptive and pose-dependent fields
    (canonical text, pose_rule, face selectors, origins, bounding box) are left out.
    """
    fr, geo = sc["frames"], sc["geometry"]
    return {
        "frames": {k: fr[k] for k in ("hfss_global_axes_in_canonical", "exit_cs_axes_in_canonical")},
        "symmetry": sc["symmetry"], "boundaries": sc["boundaries"],
        "geometry": {k: geo[k] for k in ("shape", "extent_mm") if k in geo},
        "modes": {k: v for k, v in sc["modes"].items() if k not in ("propagating_count", "basis")},
        "cross_section": sc["exit_field"]["cross_section"],
    }


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
    return [k for k in ia if not _close(ia[k], ib[k])]


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


def check(sc, dataset_id, stem, dir_freq_ghz):
    """F1-F11 and F13, as BBRDatasetSidecar::Parse runs them."""
    if not isinstance(sc, dict):
        raise ValueError("BBR024: the top level must be a JSON object")
    schema = _req(sc, "schema_version", "", "str")                                # F1
    if schema != "1" and not schema.startswith("1."):
        raise ValueError(f"BBR024: schema_version {schema}: BBRsim reads schema 1.x")
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
    if _req(xf, "columns", "exit_field.", "list") != WG_COLUMNS:
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

    tr = _req(sc, "transmittance", "", "dict")                                     # F9
    _expect(_req(tr, "definition", "transmittance.", "str"), TRANSMITTANCE, "transmittance.definition")
    if _req(tr, "incoming_includes_cos_theta", "transmittance.", "bool"):
        raise ValueError("BBR025: transmittance.incoming_includes_cos_theta must be false")

    sym = _req(sc, "symmetry", "", "dict")                                         # F10
    for k in ("mirror_l", "mirror_g", "end_to_end"):
        if not _req(sym, k, "symmetry.", "bool"):
            raise ValueError(f"BBR025: symmetry.{k} is false; the sampler folds by both mirrors "
                             "and serves both ends from one table")
    _req(sym, "rotational", "symmetry.", "bool")
    _req(sc, "boundaries", "", "dict")

    ext = _req(_req(sc, "geometry", "", "dict"), "extent_mm", "geometry.", "dict")  # F11
    for k in "plg":
        if not _req(ext, k, "geometry.extent_mm.", "num") > 0:
            raise ValueError(f"BBR024: geometry.extent_mm.{k} must be > 0")

    md = _req(sc, "modes", "", "dict")                                             # F13
    got = (_req(md, "mode", "modes.", "str"), _req(md, "cutoff_ghz", "modes.", "num"),
           _req(md, "polarization_filter_limit_ghz", "modes.", "num"))
    want = lowest_mode(cs)
    if got[0] != want[0] or abs(got[1] - want[1]) > MODE_REL * want[1] or abs(got[2] - want[2]) > MODE_REL * want[2]:
        raise ValueError(f"BBR025: modes (mode, cutoff_ghz, polarization_filter_limit_ghz) = {got}, "
                         f"the declared cross-section gives {want}")


def _contains(section, y, z):
    if section["shape"] == "rectangle":
        return (abs(y) <= section["y_e_half_m"] * (1 + SECTION_REL)
                and abs(z) <= section["z_e_half_m"] * (1 + SECTION_REL))
    return y * y + z * z <= section["radius_m"] ** 2 * (1 + SECTION_REL)


def _check_one(stem, d, sc, exempt):
    """C1-C5 on the description of one polarization; exempt holds the (key, Y, Z) that C5 lets
    lie outside the section."""
    tag = f"{stem}_Ephi={d['ephi']}"
    ff, xf = sc["far_field"], sc["exit_field"]
    if d["far_field"]["columns"] != ff["columns"] or d["exit_field"]["columns"] != xf["columns"]:
        raise ValueError(f"BBR013: {tag}: a CSV header differs from the sidecar's columns")              # C1
    # C2: keys rounded to 0.01 deg, as BBRHFSSData rounds them; the waveguide keys must lie on the
    # declared grid, and the far-field keys must equal them (BBR007 in BBRHFSSData too)
    phis = {round(float(v), 2) for v in sc["excitation"]["incident_phi_deg"]}
    thetas = {round(float(v), 2) for v in sc["excitation"]["incident_theta_deg"]}
    off = sorted(k for k in d["keys"] if k[0] not in phis or k[1] not in thetas)
    if off:
        raise ValueError(f"BBR007: {tag}: keys {off} are not on the declared incidence grid")
    if d["ff_keys"] != d["keys"]:
        raise ValueError(f"BBR007: {tag}: far_field.csv and waveguide.csv disagree on incidence keys "
                         f"(far field only {sorted(d['ff_keys'] - d['keys'])}, "
                         f"waveguide only {sorted(d['keys'] - d['ff_keys'])})")
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


def check_csvs(base, stem, sc):
    """C1-C5 on both polarizations, as BBRHFSSData runs them when it is given the sidecar, plus
    the far-field step and the agreement of Ephi=1 with Ephi=0 on keys, grids and per-key counts.

    Precondition: sc has passed check(); its blocks and types are not re-checked here. A missing
    far_field.csv or waveguide.csv raises ValueError starting BBR001 or BBR002, as in BBRHFSSData.
    """
    d0, d1 = describe_csvs(base, stem)["polarizations"]
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
    """Everything: check(), the full modes block (basis is recorded only, not compared), the
    sha256, bytes and rows of the four CSVs, and check_csvs(). Unknown fields are ignored."""
    sc = load(base, stem)
    check(sc, dataset_id, stem, dir_freq_ghz)
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
    return sc
