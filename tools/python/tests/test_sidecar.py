"""bbrsim.sidecar: mode derivations, the schema checks and the CSV checks."""
import copy
import math
import os

import pytest

from bbrsim import sidecar

FF_HDR = ",".join(sidecar.FF_COLUMNS) + "\n"
WG_HDR = ",".join(sidecar.WG_COLUMNS) + "\n"


def write_tiny(base, dataset_id="c", label="500GHz", exit_yz=(1e-3, 0.0), field=("0,1,0,0,0,0",) * 2):
    """One key (0, 180), two far-field rows, one exit point; T0 = 0.5, T1 = 0.25.
    field: the six exit-field components of the exit point, per polarization."""
    stem = f"{dataset_id}_{label}"
    for e, out in ((0, 0.5), (1, 0.25)):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d)
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(FF_HDR + "".join(f"{label},{e},0,180,0,{t},0.1,0.0,0.2,0.05\n" for t in (60, 120)))
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(WG_HDR + f"{label},{e},0,180,{out},1.0,0,{exit_yz[0]},{exit_yz[1]},{field[e]}\n")
    return stem


def tiny_sidecar(base, stem, dataset_id="c", label="500GHz", half=(2e-3, 2.5e-5)):
    return sidecar.build_from_csvs(
        base, dataset_id, label,
        section={"shape": "rectangle", "y_e_half_m": half[0], "z_e_half_m": half[1]},
        extent_mm={"p": 1.0, "l": 4.0, "g": 0.05}, provenance={"producer": "test"},
        exit_origin_mm=[0.0, 0.0, 1.0], plane_wave_origin_mm=[0.0, 0.0, 0.0],
        bounding_box_mm=[-0.025, -2.0, 0.0, 0.025, 2.0, 1.0])


def test_rect_modes_crack1():
    m = sidecar.rect_modes(0.01, 5e-05, 500.0)
    assert m["mode"] == "TE10" and m["cutoff_ghz"] == pytest.approx(14.9896229, rel=1e-9)
    assert m["polarization_filter_limit_ghz"] == pytest.approx(2997.92458, rel=1e-9)
    assert [o["mode"] for o in m["gap_family_onsets"]] == ["TE10", "TE01", "TE02", "TE03", "TE04", "TE05", "TE06"]
    assert m["mode_count_below_limit"] == 13970 and m["propagating_count"] == 33


def test_rect_modes_crack2_labels():
    m = sidecar.rect_modes(0.01, 1e-04, 500.0)
    labels = [o["mode"] for o in m["gap_family_onsets"]]
    assert len(labels) == 14 and labels[9:11] == ["TE09", "TE0,10"]
    assert m["mode_count_below_limit"] == 28061 and m["propagating_count"] == 33
    assert m["polarization_filter_limit_ghz"] == pytest.approx(1498.96229, rel=1e-9)


def test_disc_modes_round_gap():
    m = sidecar.disc_modes(5e-05, 2000.0)
    assert m["mode"] == "TE11" and m["cutoff_ghz"] == pytest.approx(1756.98, abs=0.01)
    assert len(m["cutoffs"]) == 114 and m["propagating_count"] == 1
    head = [(c["mode"], round(c["cutoff_ghz"], 2), c["degeneracy"]) for c in m["cutoffs"][:10]]
    assert head == [("TE11", 1756.98, 2), ("TM01", 2294.85, 1), ("TE21", 2914.56, 2), ("TE01", 3656.48, 1),
                    ("TM11", 3656.48, 2), ("TE31", 4009.06, 2), ("TM21", 4900.77, 2), ("TE41", 5074.38, 2),
                    ("TE12", 5087.63, 2), ("TM02", 5267.64, 1)]
    assert any("," in c["mode"] for c in m["cutoffs"])          # n reaches 10 within 20 THz


def test_mode_label():
    assert sidecar.mode_label("TE", 1, 1) == "TE11"
    assert sidecar.mode_label("TE", 18, 1) == "TE18,1" and sidecar.mode_label("TM", 0, 10) == "TM0,10"


def test_build_then_full_check(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sidecar.write(base, stem, sc)
    assert sidecar.check_full(base, "c", stem, 500.0)["modes"]["mode"] == "TE10"


@pytest.mark.parametrize("edit,code", [
    (lambda s: s.update(schema_version="2.0"), "BBR024"),
    (lambda s: s.update(dataset_id="other"), "BBR024"),
    (lambda s: s.update(frequency_label="501GHz"), "BBR024"),
    (lambda s: s["far_field"].update(points_per_key="2"), "BBR024"),
    (lambda s: s["frames"].update(hfss_global_axes_in_canonical={"x": [0, 0, 1], "y": [0, -1, 0], "z": [1, 0, 0]}), "BBR025"),
    (lambda s: s["frames"]["exit_cs"].update(z=[1, 0, 0]), "BBR025"),
    (lambda s: s["excitation"].update(incidence_convention="propagation_direction"), "BBR025"),
    (lambda s: s["exit_field"].update(field_components_frame="exit_cs"), "BBR025"),
    (lambda s: s["exit_field"].update(cross_section={"shape": "polygon", "vertices_m": [[0, 0]]}), "BBR025"),
    (lambda s: s["transmittance"].update(incoming_includes_cos_theta=True), "BBR025"),
    (lambda s: s["symmetry"].update(mirror_l=False), "BBR025"),
    (lambda s: s["modes"].update(cutoff_ghz=15.5), "BBR025"),
    (lambda s: s["far_field"].update(columns=s["far_field"]["columns"][::-1]), "BBR025"),
], ids=["schema", "id", "label", "type", "rotated-frame", "exit-z", "incidence", "field-frame",
        "polygon", "cos-theta", "mirror", "cutoff", "columns"])
def test_check_rejects(tmp_path, edit, code):
    base = str(tmp_path)
    stem = write_tiny(base)
    bad = copy.deepcopy(tiny_sidecar(base, stem))
    edit(bad)
    with pytest.raises(ValueError, match=f"^{code}"):
        sidecar.check(bad, "c", stem, 500.0)


def test_check_accepts_unknown_fields(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)
    ok = tiny_sidecar(base, stem)
    ok["x_future_field"] = {"a": 1}
    ok["schema_version"] = "1.4"
    ok["frequency_ghz"] = 500                       # an integer is a number
    sidecar.check(ok, "c", stem, 500.0)
    ok["modes"]["x_future"] = 1                     # nested unknown fields: check_full ignores them too
    ok["modes"]["gap_family_onsets"][0]["x_future"] = "a"
    for entry in ok["files"].values():
        entry["mtime_utc"] = "2026-10-05T00:00:00Z"
    ok["files"]["c_500GHz_Ephi=0/notes.txt"] = {"sha256": "0" * 64}
    sidecar.write(base, stem, ok)
    sidecar.check_full(base, "c", stem, 500.0)


def test_check_csvs_rejects_point_outside_section(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base, exit_yz=(3e-3, 0.0))     # 3 mm off axis
    sc = tiny_sidecar(base, stem)                     # section half 2 mm
    with pytest.raises(ValueError, match="^BBR025"):
        sidecar.check_csvs(base, stem, sc)


def test_check_full_rejects_stale_checksum(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sidecar.write(base, stem, sc)
    with open(os.path.join(base, f"{stem}_Ephi=0", "far_field.csv"), "a") as fh:
        fh.write("\n")
    with pytest.raises(ValueError, match="^BBR024"):
        sidecar.check_full(base, "c", stem, 500.0)


# --- fix round 1: both polarizations, forward compatibility, rows, steps, outside points ---------

def rewrite(base, stem, ephi, name, transform):
    p = os.path.join(base, f"{stem}_Ephi={ephi}", name)
    with open(p) as fh:
        text = fh.read()
    with open(p, "w") as fh:
        fh.write(transform(text))


def write_with_files(base, stem, sc):
    """Refresh the files block from the CSVs on disk, so check_full gets past the checksums."""
    sc["files"] = sidecar.file_entries(base, stem)
    sidecar.write(base, stem, sc)


def widen_y_grid(sc):
    sc["exit_field"]["grid"]["y_e"] = {"min": 0.001, "max": 0.009, "count": 1}


@pytest.mark.parametrize("name,transform,match,edit", [
    ("waveguide.csv", lambda t: "A,B,C\n" + t.split("\n", 1)[1], r"BBR013: .*Ephi=1.*header", None),
    ("waveguide.csv", lambda t: t.replace(",0,180,", ",45,135,"), r"BBR007: c_500GHz_Ephi=1: keys", None),
    ("waveguide.csv", lambda t: t.replace(",1.0,0,", ",1.0,0.5,"), r"BBR012: c_500GHz_Ephi=1: exit points off", None),
    ("waveguide.csv", lambda t: t.replace(",0.001,0.0,", ",0.009,0.0,"), r"BBR025: c_500GHz_Ephi=1: 1 exit point",
     widen_y_grid),
    ("far_field.csv", lambda t: t.replace(",0,120,", ",0,90,"), r"BBR012: c_500GHz_Ephi=1: far-field theta_deg step",
     None),
], ids=["header", "off-grid-key", "x-off-plane", "outside-section", "far-field-step"])
def test_check_full_checks_ephi1(tmp_path, name, transform, match, edit):
    """A defect in the Ephi=1 CSVs alone is caught: C1-C5 run on both polarizations."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    if edit:
        edit(sc)                                    # so that the 9 mm point passes C4 and reaches C5
    rewrite(base, stem, 1, name, transform)
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=f"^{match}"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_check_csvs_polarizations_agree(tmp_path):
    """Ephi=1 keys on the declared grid but unlike Ephi=0's: BBR012."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["excitation"]["incident_theta_deg"] = [135.0, 180.0]
    for name in ("far_field.csv", "waveguide.csv"):
        rewrite(base, stem, 1, name, lambda t: t.replace(",0,180,", ",0,135,"))
    with pytest.raises(ValueError, match="^BBR012: c_500GHz: Ephi=1 incidence keys"):
        sidecar.check_csvs(base, stem, sc)


def test_check_csvs_rejects_wrong_step(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    assert sc["far_field"]["theta_deg"]["step"] == 60.0
    sc["far_field"]["theta_deg"]["step"] = 30.0
    with pytest.raises(ValueError, match="^BBR012: .*theta_deg step 30"):
        sidecar.check_csvs(base, stem, sc)


@pytest.mark.parametrize("name,code", [("far_field.csv", "BBR001"), ("waveguide.csv", "BBR002")])
def test_check_csvs_missing_csv(tmp_path, name, code):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    os.remove(os.path.join(base, f"{stem}_Ephi=1", name))
    with pytest.raises(ValueError, match=f"^{code}: .*{name}"):
        sidecar.check_csvs(base, stem, sc)


def test_check_full_ignores_basis(tmp_path):
    """modes.basis is recorded only: the HFSS writer words it its own way."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["modes"]["basis"] = "closed ideal-PEC rectangular guide, worded by another writer"
    sidecar.write(base, stem, sc)
    sidecar.check_full(base, "c", stem, 500.0)


def test_file_entries_rows(tmp_path):
    """rows = non-blank lines after the header, as pandas counts them."""
    base = str(tmp_path)
    stem = write_tiny(base)
    rewrite(base, stem, 0, "waveguide.csv", lambda t: t.rstrip("\n"))              # no trailing newline
    rewrite(base, stem, 1, "waveguide.csv", lambda t: t + "\n  \n")               # trailing blank lines
    rewrite(base, stem, 0, "far_field.csv", lambda t: t.replace("\n", "\r\n"))    # CRLF
    rows = {k: v["rows"] for k, v in sidecar.file_entries(base, stem).items()}
    assert rows == {f"{stem}_Ephi=0/far_field.csv": 2, f"{stem}_Ephi=0/waveguide.csv": 1,
                    f"{stem}_Ephi=1/far_field.csv": 2, f"{stem}_Ephi=1/waveguide.csv": 1}


@pytest.mark.parametrize("label,dir_freq,match", [
    ("xGHz", 500.0, "frequency_label 'xGHz' does not parse"),
    ("501GHz", 500.0, "frequency_label '501GHz' .* disagrees"),
], ids=["unparsable", "label-value"])
def test_check_f3_names_the_label(tmp_path, label, dir_freq, match):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["frequency_label"] = label
    with pytest.raises(ValueError, match=f"^BBR024: {match}"):
        sidecar.check(sc, "c", f"c_{label}", dir_freq)


def test_check_f3_names_frequency_ghz(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["frequency_ghz"] = 510.0
    with pytest.raises(ValueError, match="^BBR024: frequency_ghz 510.0 disagrees"):
        sidecar.check(sc, "c", stem, 500.0)


ZERO = "0,0,0,0,0,0"


@pytest.mark.parametrize("outside_points,field,ok", [
    ("zero", (ZERO, ZERO), True),                     # zero in both polarizations: exempt
    ("zero", (ZERO, "0,1,0,0,0,0"), False),           # nonzero in Ephi=1: not exempt
    ("none", (ZERO, ZERO), False),                    # every point must lie inside
    ("omitted", (ZERO, ZERO), False),
], ids=["zero-exempt", "zero-nonzero-ephi1", "none", "omitted"])
def test_check_csvs_outside_points(tmp_path, outside_points, field, ok):
    base = str(tmp_path)
    stem = write_tiny(base, exit_yz=(3e-3, 0.0), field=field)    # 3 mm off axis, section half 2 mm
    sc = tiny_sidecar(base, stem)
    sc["exit_field"]["outside_points"] = outside_points
    sidecar.check(sc, "c", stem, 500.0)
    if ok:
        sidecar.check_csvs(base, stem, sc)
    else:
        with pytest.raises(ValueError, match="^BBR025: .*outside the declared cross-section"):
            sidecar.check_csvs(base, stem, sc)


# --- fix round 2: the physics-only invariant, the disc modes, load(), malformed CSVs -------------

def built(tmp_path):
    base = str(tmp_path)
    return tiny_sidecar(base, write_tiny(base))


def test_invariant_keys(tmp_path):
    assert set(sidecar.invariant(built(tmp_path))) == {"frames", "symmetry", "boundaries", "geometry",
                                                      "modes", "cross_section"}


@pytest.mark.parametrize("edit", [
    lambda s: s["frames"].update(canonical="p, l, g = propagation, long, gap (another writer's words)"),
    lambda s: s["frames"].update(pose_rule="legacy"),
    lambda s: s["frames"]["exit_cs"].update(origin_mm_global=[0.0, -5.0, 1.0]),
    lambda s: s["excitation"].update(origin_mm_global=[0.0, -5.0, 0.0]),
    lambda s: s["geometry"].update(bounding_box_mm=[-0.025, -10.0, 0.0, 0.025, 0.0, 1.0]),
    lambda s: s["modes"].update(basis="closed ideal-PEC rectangular guide, worded by another writer"),
    lambda s: s["modes"].update(propagating_count=1),
    lambda s: s.update(provenance={"producer": "hand", "hand_written": True}),
    lambda s: s.update(files={}),
    lambda s: s["symmetry"].update(x_future=True),
], ids=["canonical", "pose-rule", "exit-origin", "plane-wave-origin", "bounding-box", "basis",
        "propagating-count", "provenance", "files", "extra-symmetry-key"])
def test_same_invariant_ignores_descriptive_fields(tmp_path, edit):
    a = built(tmp_path)
    b = copy.deepcopy(a)
    edit(b)
    assert sidecar.same_invariant(a, b) and sidecar.same_invariant(b, a)


@pytest.mark.parametrize("edit", [
    lambda s: s["geometry"]["extent_mm"].update(p=4.0),
    lambda s: s["exit_field"]["cross_section"].update(y_e_half_m=2.5e-3),
    lambda s: s["symmetry"].update(rotational=True),
    lambda s: s["boundaries"].update(walls="finite-conductivity"),
    lambda s: s["modes"].update(cutoff_ghz=s["modes"]["cutoff_ghz"] * (1 + 1e-6)),
    lambda s: s["modes"].update(gap_family_onsets=s["modes"]["gap_family_onsets"][:-1]),
], ids=["extent", "cross-section", "symmetry-flag", "walls", "cutoff-1e-6", "onsets-shorter"])
def test_same_invariant_detects_physics(tmp_path, edit):
    a = built(tmp_path)
    b = copy.deepcopy(a)
    edit(b)
    assert not sidecar.same_invariant(a, b) and not sidecar.same_invariant(b, a)


def test_same_invariant_numbers_to_1e9(tmp_path):
    """The two writers' cutoffs differ by about 2e-11 relative: the same."""
    a = built(tmp_path)
    b = copy.deepcopy(a)
    b["modes"]["cutoff_ghz"] *= 1 + 1e-12
    assert sidecar.same_invariant(a, b)


def test_modes_for_disc():
    r = 5e-05
    m = sidecar.modes_for({"shape": "disc", "radius_m": r}, 2000.0)
    cut = sidecar.X11_PRIME * sidecar.C / (2 * math.pi * r) / 1e9
    assert m["mode"] == "TE11" and m["cutoff_ghz"] == pytest.approx(cut, rel=1e-12)
    assert m["polarization_filter_limit_ghz"] == m["cutoff_ghz"]
    assert m["basis"].startswith("closed ideal-PEC circular guide, R = 5e-05 m")
    mode, low, filt = sidecar.lowest_mode({"shape": "disc", "radius_m": r})
    assert mode == m["mode"] and low == pytest.approx(m["cutoff_ghz"], rel=1e-12)
    assert filt == pytest.approx(m["polarization_filter_limit_ghz"], rel=1e-12)


def test_load_accepts_bom(tmp_path):
    base = str(tmp_path)
    with open(sidecar.path_for(base, "c_500GHz"), "wb") as fh:
        fh.write(b'\xef\xbb\xbf{"schema_version": "1.0"}\r\n')
    assert sidecar.load(base, "c_500GHz") == {"schema_version": "1.0"}


def test_load_missing(tmp_path):
    with pytest.raises(ValueError, match=r"^BBR024: .*c_500GHz\.dataset\.json: missing"):
        sidecar.load(str(tmp_path), "c_500GHz")


def test_load_invalid_json(tmp_path):
    base = str(tmp_path)
    with open(sidecar.path_for(base, "c_500GHz"), "w") as fh:
        fh.write('{"schema_version": ')
    with pytest.raises(ValueError, match=r"^BBR024: .*: not valid JSON"):
        sidecar.load(base, "c_500GHz")


def test_check_csvs_far_field_key_off_grid(tmp_path):
    """A key in far_field.csv only (two rows, so the per-key count holds): BBR007, as BBRHFSSData."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 0, "far_field.csv",
            lambda t: t + "".join(f"500GHz,0,45,135,0,{th},0.1,0.0,0.2,0.05\n" for th in (60, 120)))
    with pytest.raises(ValueError, match=r"^BBR007: c_500GHz_Ephi=0: far_field.csv and waveguide.csv "
                                         r"disagree .*far field only \[\(45.0, 135.0\)\]"):
        sidecar.check_csvs(base, stem, sc)


@pytest.mark.parametrize("transform,match", [
    (lambda t: t + "500GHz,1,0,180,0,90,0.1,0.0,0.2,0.05,7\n", "wrong number of fields"),   # 2nd row long
    (lambda t: t.replace(",0.05\n", ",0.05,7\n"), "more fields than the header"),            # every row long
    (lambda t: t + "500GHz,1,0,180,0\n", "lacks fields"),                                     # short row
], ids=["one-row-long", "every-row-long", "short-row"])
def test_check_csvs_malformed_row(tmp_path, transform, match):
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 1, "far_field.csv", transform)
    with pytest.raises(ValueError, match=f"^BBR013: .*Ephi=1.*far_field.csv: .*{match}"):
        sidecar.check_csvs(base, stem, sc)
