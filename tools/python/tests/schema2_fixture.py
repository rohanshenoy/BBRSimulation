"""A small schema-2.0 dataset for the tests, the Python twin of tests/HFSSFixture.hh (NineKeyGrid2,
Schema2Params, SidecarJson2): the same keys, lattice, powers, field values and sidecar content, written
here independently of bbrsim.sidecar so that the checks and the builder are tested against a reference
they do not share code with (the closed forms below are copies, as the C++ fixture's are).

Keys: IWaveTheta {90, 135, 180} x IWavePhi {0, 45, 90}. Exit lattice y {-1e-3, 0, 1e-3} m x
z {-2.5e-5, 0, 2.5e-5} m (the disc variant: +-2e-5 m on both axes, inside R = 50 um). IngoingPower
1e-9 W, t0 = 0.6, t1 = 0.3 on every key except Ephi=1 at (0, 180), which is dark (OutgoingPower 0, the
stored fields not). The six field values of exit point i are (i+1 ... i+6) * 1e-3 for Ephi=0 and
(i+6 ... i+1) * 1e-3 for Ephi=1, so every column of the 2.0 position map differs. InterceptFactor is the
key's c_k. The far field is a complete Phi {-40, 10, 60} x Theta {40, 75, 110} grid per key.
"""
import math
import os

STEM = "crack_500GHz"
FF_COLUMNS = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta",
              "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"]
WG_COLUMNS_2 = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower", "X", "Y", "Z",
                "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor"]
FF_HDR = ",".join(FF_COLUMNS) + "\n"
WG_HDR2 = ",".join(WG_COLUMNS_2) + "\n"
PHIS = (0, 45, 90)
THETAS = (90, 135, 180)
C = 299792458.0

RECT_Y, RECT_Z = (-1e-3, 0.0, 1e-3), (-2.5e-5, 0.0, 2.5e-5)
DISC_Y, DISC_Z = (-2e-5, 0.0, 2e-5), (-2e-5, 0.0, 2e-5)
FF_THETAS, FF_PHIS = ("40", "75", "110"), ("-40", "10", "60")


# --- the HF-049 rule, a fixture-local copy of the closed form -------------------------------------

def _delta(thetas):
    has90 = any(abs(t - 90.0) <= 1e-9 for t in thetas)
    above = [t for t in thetas if abs(t - 90.0) > 1e-9 and t > 90.0]
    return min(above) - 90.0 if has90 and above else -1.0


def cell90(thetas):
    d = _delta(thetas)
    if d < 0:
        return -1.0
    a0 = math.pi / 2 - d / 2 * math.pi / 180.0
    af = math.acos(0.01)
    prim = lambda a: math.log(1.0 / math.cos(a) + math.tan(a))
    length = math.pi / 2 - a0
    inv = (prim(af) - prim(a0)) + (math.pi / 2 - af) / 0.01 if af > a0 else (math.pi / 2 - a0) / 0.01
    return length / inv


def intercept(theta, thetas):
    return cell90(thetas) if abs(theta - 90.0) <= 1e-9 else abs(math.cos(theta * math.pi / 180.0))


# --- CSVs -----------------------------------------------------------------------------------------

def _fmt(v):
    return repr(float(v)) if not isinstance(v, int) else str(v)


def far_field_text(e, label="500GHz", phis=PHIS, thetas=THETAS):
    """The nine-key-grid far field of NineKeyGrid, for polarization e."""
    out = FF_HDR
    for ph in phis:
        for t in thetas:
            for th in FF_THETAS:
                for pm in FF_PHIS:
                    pt = f"{pm},{th}"
                    a0, a1 = "0.05,0.02,0.03,-0.01", "0.02,-0.03,0.04,0.01"
                    if pt == "10,40":
                        a0, a1 = "0.9,0.0,0.1,0.4", "0.3,0.1,0.8,-0.2"
                    if pt == "-40,75":
                        a0, a1 = "0.2,-0.3,0.7,0.1", "-0.5,0.2,0.1,0.6"
                    if pt == "60,110":
                        a0, a1 = "0.5,0.1,-0.6,0.3", "0.4,0.4,-0.3,0.2"
                    out += f"{label},{e},{ph},{t},{pt},{a1 if e else a0}\n"
    return out


def waveguide_rows(e, ys, zs, thetas=THETAS, phis=PHIS, label="500GHz", dark_key=(0, 180)):
    """The 16-column rows of polarization e. tp = 0.6 (Ephi=0) or 0.3 (Ephi=1) of A.I = 1e-9 W."""
    out = WG_HDR2
    for ph in phis:
        for t in thetas:
            dark = e == 1 and (ph, t) == dark_key
            tp = 0.6 if e == 0 else 0.3
            ck = repr(intercept(t, thetas))
            i = 0
            for y in ys:
                for z in zs:
                    fields = ",".join(_fmt((i + 1 + f if e == 0 else i + 6 - f) * 1e-3) for f in range(6))
                    out += f"{label},{e},{ph},{t},{_fmt(0.0 if dark else tp * 1e-9)},1e-09,0,{_fmt(y)},{_fmt(z)},{fields},{ck}\n"
                    i += 1
    return out


def write_nine_key_grid2(base, stem=STEM, disc=False, thetas=THETAS, phis=PHIS):
    """Write <stem>_Ephi={0,1}/{far_field,waveguide}.csv under base; return stem."""
    ys, zs = (DISC_Y, DISC_Z) if disc else (RECT_Y, RECT_Z)
    label = stem.rsplit("_", 1)[1]
    for e in (0, 1):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(far_field_text(e, label, phis, thetas))
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(waveguide_rows(e, ys, zs, thetas, phis, label))
    return stem


# --- the sidecar of that dataset ------------------------------------------------------------------

def sidecar2(stem=STEM, *, disc=False, radius=5e-5, y_half=4.5e-3, z_half=2.5e-5, thetas=THETAS, phis=PHIS,
             tflux=None, emitter="uniform_theta", floor=0.01):
    """The 2.0 sidecar of write_nine_key_grid2, as a dict (SidecarJson2 of the C++ fixture). tflux: one
    [t0, t1, dark0, dark1] per key in phi-major, then theta order; default 0.6, 0.3 and (0, 180)
    Ephi=1 dark."""
    label = stem.rsplit("_", 1)[1]
    ys, zs = (DISC_Y, DISC_Z) if disc else (RECT_Y, RECT_Z)
    if tflux is None:
        tflux = [[0.6, 0.0 if (ph, t) == (0, 180) else 0.3, 0, 1 if (ph, t) == (0, 180) else 0]
                 for ph in phis for t in thetas]
    section = ({"shape": "disc", "radius_m": radius} if disc else
               {"shape": "rectangle", "y_e_half_m": y_half, "z_e_half_m": z_half})
    if disc:
        cut = 1.8411837813 * C / (2 * math.pi * radius) / 1e9
        mode, filt = "TE11", cut
        regime = [{"mode": n, "cutoff_ghz": x * C / (2 * math.pi * radius) / 1e9}
                  for n, x in (("TE11", 1.8411837813), ("TM01", 2.404825557695773), ("TE21", 3.054236928227140))]
    else:
        a, b = 2 * y_half, 2 * z_half
        mode, cut = ("TE10", C / (2 * a) / 1e9) if a >= b else ("TE01", C / (2 * b) / 1e9)
        filt = C / (2 * b) / 1e9
        regime = [{"mode": "TE10", "cutoff_ghz": C / (2 * a) / 1e9}]
        n = 1
        while n * C / (2 * b) / 1e9 <= 20000.0 + 1e-9:
            regime.append({"mode": f"TE0{n}" if n < 10 else f"TE0,{n}", "cutoff_ghz": n * C / (2 * b) / 1e9})
            n += 1
        regime.sort(key=lambda r: r["cutoff_ghz"])
    delta = _delta(thetas)
    c_key, per_key, scale, ff = [], [], [], []
    k = 0
    for ph in phis:
        for t in thetas:
            t0, t1, d0, d1 = tflux[k]
            ck = intercept(t, thetas)
            c_key.append({"phi_deg": ph, "theta_deg": t, "c_k": ck})
            per_key.append({"phi_deg": ph, "theta_deg": t, "c_k": ck, "polarizations": {
                "0": {"t_flux": 0.0 if d0 else t0, "outgoing_power_abs_re_s_w": 1.05e-9 * t0},
                "1": {"t_flux": 0.0 if d1 else t1, "outgoing_power_abs_re_s_w": 1.05e-9 * t1}}})
            scale.append({"phi_deg": ph, "theta_deg": t, "w0": 1.0 if d0 else 0.9, "w1": 1.0 if d1 else 0.8,
                          "dark0": bool(d0), "dark1": bool(d1)})
            ff.append({"phi_deg": ph, "theta_deg": t,
                       "polarizations": {"0": None if d0 else 0.179, "1": None if d1 else 0.179}})
            k += 1
    ny, nz = len(ys), len(zs)
    # extent_mm l and g: the section's full widths (a disc: its diameter), as the C++ SidecarJson
    ext_l, ext_g = (2e3 * radius, 2e3 * radius) if disc else (2e3 * y_half, 2e3 * z_half)
    return {
        "schema_version": "2.0",
        "dataset_id": stem.rsplit("_", 1)[0],
        "frequency_ghz": float(label[:-3]),
        "frequency_label": label,
        "provenance": {"producer": "tests/schema2_fixture.py"},
        "model": {"kind": "screen_free_abc"},
        "conventions": {"time_dependence": "exp(+j omega t)", "phasor": "peak"},
        "frames": {
            "canonical": "p,l,g", "pose_rule": "canonical-z",
            "hfss_global_axes_in_canonical": {"x": [0, 0, -1], "y": [0, 1, 0], "z": [1, 0, 0]},
            "entrance_face": {"axis": "z", "side": "min"}, "exit_face": {"axis": "z", "side": "max"},
            "entrance_outward_normal_global": [0, 0, -1], "exit_outward_normal_global": [0, 0, 1],
            "exit_cs": {"name": "outgoing_cs", "origin": "exit_face_center", "origin_mm_global": [0, 0, 4],
                        "x": [0, 0, 1], "y": [0, 1, 0], "z": [-1, 0, 0]},
            "exit_cs_axes_in_canonical": {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]},
        },
        "excitation": {"coordinate_system": "global", "incidence_convention": "arrival_direction",
                       "normal_entry_theta_deg": 180, "incident_phi_deg": list(phis),
                       "incident_theta_deg": list(thetas),
                       "polarization_convention": "Ephi=0: E_theta=1; Ephi=1: E_phi=1",
                       "origin_mm_global": [0, 0, 0]},
        "far_field": {"coordinate_system": "exit_cs", "definition": "Theta-Phi",
                      "component_basis": "spherical_in_exit_cs",
                      "theta_deg": {"min": 40.0, "max": 110.0, "count": 3, "step": 35.0},
                      "phi_deg": {"min": -40.0, "max": 60.0, "count": 3, "step": 50.0},
                      "points_per_key": 9, "columns": list(FF_COLUMNS),
                      "normalization": "pattern_only", "ff_over_flux": ff},
        "exit_field": {
            "coordinate_system": "exit_cs", "points_in_si": True, "field_in_ref_cs": False,
            "field_components_frame": "hfss_global", "plane": "x_e=0",
            "grid": {"y_e": {"min": ys[0], "max": ys[-1], "count": ny, "step": (ys[-1] - ys[0]) / (ny - 1)},
                     "z_e": {"min": zs[0], "max": zs[-1], "count": nz, "step": (zs[-1] - zs[0]) / (nz - 1)}},
            "cross_section": section, "outside_points": "none", "rim_points": "included",
            "points_per_key_retained": ny * nz, "columns": list(WG_COLUMNS_2),
            "tangential_components": ["Ex", "Ey"],
            "quadrature": "uniform_disc" if disc else "trapezoid",
            "sampler": "abs_e_t_squared_times_quadrature",
            "scale": {"rule": "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w", "per_key": scale},
            "column_order_check": {"rule": "reference_exporter_interleaved_v1", "method": "reference",
                                   "reference": {"dataset_id": "InfParallelPlate_crack1Rohan",
                                                 "frequency_ghz": 500, "passed": True}},
            "plane_tolerance_m": 1e-12,
        },
        "transmittance": {
            "definition": "OutgoingPower / (IngoingPower * InterceptFactor)",
            "numerator": "exit_flux_re_s_dot_n", "incoming": "area_times_intensity_undivided",
            "above_one": "per_photon_cap", "c_key": c_key, "cell_90": cell90(thetas),
            "cell_90_assumes": {"emitter": emitter, "cos_floor": floor, "band_deg": [90.0 - delta / 2, 90.0]},
            "per_key": per_key},
        "cross_term": {"density": "e_t_quadrature",
                       "definition": "Re rho = sum_i w_i Re(E_t0,i . conj E_t1,i) / sqrt(sum_i w_i abs(E_t0,i)^2 "
                                     "* sum_i w_i abs(E_t1,i)^2)",
                       "valid_below_ghz": filt, "above_limit": "warn_once"},
        "symmetry": {"mirror_l": True, "mirror_g": True, "end_to_end": True, "rotational": bool(disc)},
        "boundaries": {"entrance": "radiation", "exit": "radiation", "walls": "PEC"},
        "geometry": {"shape": "cylinder" if disc else "box", "extent_mm": {"p": 4, "l": ext_l, "g": ext_g},
                     "bounding_box_mm": [-0.026, -5, 0, 0.026, 5, 4]},
        "modes": {"mode": mode, "cutoff_ghz": cut, "polarization_filter_limit_ghz": filt,
                  "regime_cutoffs_ghz": regime},
    }
