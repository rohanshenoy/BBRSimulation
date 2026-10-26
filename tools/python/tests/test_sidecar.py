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


@pytest.mark.parametrize("radius_m,modes", [
    (4.39e-6, []),                     # TE11 just above the 20 THz list limit (at R = 4.39246 um it is on it)
    (4.40e-6, ["TE11"]),               # TE11 below the limit, order 0 still empty: the CDX-006 window
    (5.73e-6, ["TE11"]),               # TM01, the first order-0 root, just above the limit
    (5.75e-6, ["TE11", "TM01"]),       # TM01 enters (at R = 5.73713 um it is on the limit)
])
def test_disc_modes_order_zero_empty(radius_m, modes):
    """CDX-006: between R = 4.39246 and 5.73713 um order 0 has no cutoff below 20 THz while TE11
    has; the search must not stop at the empty order 0."""
    assert [c["mode"] for c in sidecar.disc_modes(radius_m, 18000.0)["cutoffs"]] == modes


def test_disc_modes_5um_matches_the_producer():
    """The external review's probe: R = 5 um at 18000 GHz; the HFSS producer lists TE11 at 17569.85 GHz, count 1."""
    m = sidecar.disc_modes(5e-6, 18000.0)
    assert [(c["mode"], round(c["cutoff_ghz"], 2), c["degeneracy"]) for c in m["cutoffs"]] == [("TE11", 17569.85, 2)]
    assert m["propagating_count"] == 1


@pytest.mark.parametrize("freq_ghz,count", [(2000.0, 1), (19000.0, None), (20000.0, 114), (25000.0, 179),
                                             (30000.0, 255)])
def test_disc_propagating_count_beyond_the_list(freq_ghz, count):
    """CDX-007: the list stops at 20 THz, the count runs to the file's frequency (the HFSS producer's
    179 at 25 THz and 255 at 30 THz for R = 50 um); below the limit it counts the list."""
    m = sidecar.disc_modes(5e-5, freq_ghz)
    assert len(m["cutoffs"]) == 114
    want = count if count is not None else sum(1 for c in m["cutoffs"] if c["cutoff_ghz"] <= freq_ghz)
    assert m["propagating_count"] == want


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
    (lambda s: s.update(schema_version="3.0"), "BBR024"),
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
    (lambda s: s["boundaries"].update(walls="finite-conductivity"), "BBR025"),
    (lambda s: s["boundaries"].update(entrance="port"), "BBR025"),
    (lambda s: s["boundaries"].pop("exit"), "BBR024"),
], ids=["schema", "id", "label", "type", "rotated-frame", "exit-z", "incidence", "field-frame",
        "polygon", "cos-theta", "mirror", "cutoff", "columns", "walls", "entrance", "exit-missing"])
def test_check_rejects(tmp_path, edit, code):
    base = str(tmp_path)
    stem = write_tiny(base)
    bad = copy.deepcopy(tiny_sidecar(base, stem))
    edit(bad)
    with pytest.raises(ValueError, match=f"^{code}"):
        sidecar.check(bad, "c", stem, 500.0)

@pytest.mark.parametrize("key,factor,ok", [
    ("l", 1 + 2.5e-3, False), ("l", 1 + 2e-6, False), ("l", 1 + 5e-7, True), ("l", 1 - 5e-7, True),
    ("g", 1 + 2e-6, False), ("g", 1 - 2e-6, False), ("g", 1 + 5e-7, True),
], ids=["l-big", "l-2e-6", "l-5e-7", "l-minus-5e-7", "g-2e-6", "g-minus-2e-6", "g-5e-7"])
def test_extent_section(tmp_path, key, factor, ok):
    """CDX-028: a rectangle's extent_mm.l and g equal its declared section's full widths to 1e-6 relative."""
    base = str(tmp_path)
    stem = write_tiny(base)
    bad = copy.deepcopy(tiny_sidecar(base, stem))
    bad["geometry"]["extent_mm"][key] *= factor
    if ok:
        sidecar.check(bad, "c", stem, 500.0)
        return
    with pytest.raises(ValueError, match=rf"^BBR024: geometry\.extent_mm\.{key} .*exit_field\.cross_section"):
        sidecar.check(bad, "c", stem, 500.0)


@pytest.mark.parametrize("value,code", [
    ("vacuum", None), ("Vacuum", None), ("VACUUM", None), (None, None),
    ("silicon", "BBR025"), ("", "BBR025"),
    (3, "BBR024"), (True, "BBR024"), (["vacuum"], "BBR024"),
])
def test_material_vacuum(tmp_path, value, code):
    """CDX-029: provenance.material, when a string, is 'vacuum' (case-insensitive); not a string or null is BBR024."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = copy.deepcopy(tiny_sidecar(base, stem))
    sc["provenance"]["material"] = value
    if code is None:
        sidecar.check(sc, "c", stem, 500.0)
        return
    with pytest.raises(ValueError, match=rf"^{code}: provenance\.material") as e:
        sidecar.check(sc, "c", stem, 500.0)
    if isinstance(value, str):
        assert repr(value) in str(e.value)
    elif code == "BBR024":
        assert "must be a string or null" in str(e.value)     # schema 1 accepts null, so the message may name it


def test_material_absent_or_null_schema1(tmp_path):
    """The legacy sidecars carry material: null; null, an absent key and an absent provenance load under 1.x."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = copy.deepcopy(tiny_sidecar(base, stem))
    sc["provenance"].pop("material", None)
    sidecar.check(sc, "c", stem, 500.0)
    sc["provenance"]["material"] = None
    sidecar.check(sc, "c", stem, 500.0)
    del sc["provenance"]
    sidecar.check(sc, "c", stem, 500.0)


def test_s2_material_required():
    """CDX-029, OR-039: under 2.x provenance.material is required; null, an absent key or an absent
    provenance is BBR024 naming the field; a string follows the shared rule."""
    sc = fx.sidecar2()
    assert sc["provenance"]["material"] == "vacuum"     # the 2.0 fixture writes it
    check2(mut2(lambda s: s["provenance"].update(material="Vacuum")))
    for edit in (lambda s: s["provenance"].update(material=None),
                 lambda s: s["provenance"].pop("material"),
                 lambda s: s.pop("provenance")):
        with pytest.raises(ValueError, match=r"^BBR024: provenance\.material"):
            check2(mut2(edit))
    with pytest.raises(ValueError, match=r"^BBR025: provenance\.material"):
        check2(mut2(lambda s: s["provenance"].update(material="silicon")))
    with pytest.raises(ValueError, match=r"^BBR024: provenance\.material must be a string under schema 2\.x") as e:
        check2(mut2(lambda s: s["provenance"].update(material=3)))
    assert "or null" not in str(e.value)                  # null is refused under 2.x, so the message must not offer it


def test_s2_builder_writes_material(tmp_path):
    """Every 2.x builder caller here passes material "vacuum", so the built sidecar passes the 2.x rule."""
    base, stem = s2_dataset(tmp_path)
    assert build2(base, stem)["provenance"]["material"] == "vacuum"


def test_extent_section_disc_not_compared(tmp_path):
    """A disc's extent_mm is the HFSS bounding box; it is not compared with 2R here."""
    base = str(tmp_path)
    stem = write_disc(base)
    sc = copy.deepcopy(disc_sidecar(base, stem, "omitted"))
    sc["geometry"]["extent_mm"].update(l=0.0123, g=0.0123)
    sidecar.check(sc, "c", stem, 500.0)


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
    """Ephi=1 keys on the declared grid but unlike Ephi=0's: since C2 requires every declared key
    (CDX-021), each polarization lacks one, and Ephi=0 is checked first (BBR007)."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["excitation"]["incident_theta_deg"] = [135.0, 180.0]
    for name in ("far_field.csv", "waveguide.csv"):
        rewrite(base, stem, 1, name, lambda t: t.replace(",0,180,", ",0,135,"))
    with pytest.raises(ValueError, match=r"^BBR007: c_500GHz_Ephi=0: declared incidence keys \[\(0.0, 135.0\)\]"):
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


def write_disc(base, dataset_id="c", label="500GHz"):
    """One key (0, 180), two far-field rows, three exit points at Y = -10, 0, +10 um on Z = 0."""
    stem = f"{dataset_id}_{label}"
    for e, out, field in ((0, 0.5, "0,1,0,0,0,0"), (1, 0.25, "0,0,1,0,0,0")):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d)
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(FF_HDR + "".join(f"{label},{e},0,180,0,{t},0.1,0.0,0.2,0.05\n" for t in (60, 120)))
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(WG_HDR + "".join(f"{label},{e},0,180,{out},1.0,0,{y},0.0,{field}\n"
                                      for y in ("-1e-05", "0", "1e-05")))
    return stem


def disc_sidecar(base, stem, outside_points, radius_m=1.5e-5):
    return sidecar.build_from_csvs(
        base, "c", "500GHz", section={"shape": "disc", "radius_m": radius_m},
        extent_mm={"p": 0.4, "l": 0.03, "g": 0.03}, provenance={"producer": "test"},
        exit_origin_mm=[0.0, 0.0, 0.4], plane_wave_origin_mm=[0.0, 0.0, 0.0],
        bounding_box_mm=[-0.015, -0.015, 0.0, 0.015, 0.015, 0.4],
        resolution_mm=[0.0, 0.01, 0.01], outside_points=outside_points, bounds_method="edge_samples",
        rotational=True)


@pytest.mark.parametrize("outside_points,ok", [("omitted", True), ("none", False), ("zero", False)],
                         ids=["omitted", "none", "zero"])
def test_check_csvs_omitted_allows_a_wider_lattice(tmp_path, outside_points, ok):
    """The producer declares the whole export lattice and, with "omitted", drops the points outside
    the section, so a round face whose radius is not a step multiple leaves the outermost lattice
    columns empty: R = 15 um at a 10 um step gives a lattice Y in [-20, 20] um with 5 columns, of
    which 3 hold a retained point. The retained distinct count may then be smaller than the declared
    one, never larger; "none" and "zero" write every lattice point and keep the equality."""
    base = str(tmp_path)
    stem = write_disc(base)
    sc = disc_sidecar(base, stem, outside_points)
    assert sc["exit_field"]["grid"]["y_e"] == {"min": -1e-5, "max": 1e-5, "count": 3}
    sc["exit_field"]["grid"]["y_e"] = {"min": -2e-5, "max": 2e-5, "count": 5}
    sidecar.check(sc, "c", stem, 500.0)
    if ok:
        sidecar.check_csvs(base, stem, sc)
    else:
        with pytest.raises(ValueError, match="^BBR012: c_500GHz_Ephi=0: exit grid y_e .*unlike the declared"):
            sidecar.check_csvs(base, stem, sc)


def test_check_csvs_omitted_refuses_more_values_than_declared(tmp_path):
    base = str(tmp_path)
    stem = write_disc(base)
    sc = disc_sidecar(base, stem, "omitted")
    sc["exit_field"]["grid"]["y_e"] = {"min": -2e-5, "max": 2e-5, "count": 2}    # two columns at +-20 um
    sidecar.check(sc, "c", stem, 500.0)
    with pytest.raises(ValueError, match="^BBR012: c_500GHz_Ephi=0: exit grid y_e .*more distinct values"):
        sidecar.check_csvs(base, stem, sc)


# --- CDX-021, CDX-022: complete declared grids, per key ------------------------------------------

GRID_KEYS = ((0.0, 135.0), (0.0, 180.0), (45.0, 135.0), (45.0, 180.0))


def write_grid(base, ff=((0, 60), (0, 120)), exit_yz=((0.001, 0.0),), skip_key=None, label="500GHz"):
    """The grid fixtures of the external review's deep-debug review (CDX-021, CDX-022): keys {0, 45} x {135, 180},
    T0 = 0.75 at (0, 180) and 0.25 elsewhere, T1 = 0.25; every key holds the far-field rows ff and the
    exit rows exit_yz, in that order; skip_key leaves one key out of all four CSVs."""
    stem = f"c_{label}"
    keys = [k for k in GRID_KEYS if k != skip_key]
    for e, field in ((0, "0,1,0,0,0,0"), (1, "0,0,1,0,0,0")):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d)
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(FF_HDR + "".join(f"{label},{e},{p},{t},{fp},{ft},0.1,0.0,0.2,0.05\n"
                                      for p, t in keys for fp, ft in ff))
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(WG_HDR + "".join(
                f"{label},{e},{p},{t},{0.75 if (e, p, t) == (0, 0.0, 180.0) else 0.25},1.0,0,{y},{z},{field}\n"
                for p, t in keys for y, z in exit_yz))
    return stem


def drop_key(text, key):
    """text without the rows of incidence key (phi, theta), as write_grid formats it."""
    return "".join(line for line in text.splitlines(keepends=True) if f",{key[0]},{key[1]}," not in line)


def test_grid_complete_passes(tmp_path):
    base = str(tmp_path)
    stem = write_grid(base)
    sidecar.write(base, stem, tiny_sidecar(base, stem))
    sidecar.check_full(base, "c", stem, 500.0)


def test_c2_declared_key_missing(tmp_path):
    """CDX-021: the builder declares {0, 45} x {135, 180} from the other keys, but (0, 180) has no
    rows in any CSV; BBRsim would serve (0, 135) at normal incidence, T0 0.25 instead of 0.75."""
    base = str(tmp_path)
    stem = write_grid(base, skip_key=(0.0, 180.0))
    sc = tiny_sidecar(base, stem)
    assert sc["excitation"]["incident_phi_deg"] == [0.0, 45.0]
    assert sc["excitation"]["incident_theta_deg"] == [135.0, 180.0]
    sidecar.write(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR007: c_500GHz_Ephi=0: declared incidence keys \[\(0.0, 180.0\)\] "
                                         r"have no rows"):
        sidecar.check_full(base, "c", stem, 500.0)


@pytest.mark.parametrize("ephi,names,match", [
    (0, ("far_field.csv",), r"^BBR007: c_500GHz_Ephi=0: far_field.csv and waveguide.csv disagree"),
    (0, ("waveguide.csv",), r"^BBR007: c_500GHz_Ephi=0: far_field.csv and waveguide.csv disagree"),
    (1, ("far_field.csv", "waveguide.csv"), r"^BBR007: c_500GHz_Ephi=1: declared incidence keys \[\(0.0, 180.0\)\]"),
], ids=["one-table-far-field", "one-table-waveguide", "one-polarization"])
def test_c2_key_missing_from_part_of_the_data(tmp_path, ephi, names, match):
    """CDX-021's other two cases: a key missing from one table, and from one polarization only."""
    base = str(tmp_path)
    stem = write_grid(base)
    sc = tiny_sidecar(base, stem)
    for name in names:
        rewrite(base, stem, ephi, name, lambda t: drop_key(t, (0.0, 180.0)))
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=match):
        sidecar.check_full(base, "c", stem, 500.0)


@pytest.mark.parametrize("kw,match", [
    (dict(ff=((0, 60), (0, 60), (30, 60), (30, 120))),
     r"far-field point \(Phi, Theta\) = \(0.0, 60.0\) appears twice in key \(0.0, 135.0\)"),
    (dict(ff=((0, 60), (0, 120), (0, 61), (30, 60), (30, 120), (30, 61))),
     r"far-field row \(Phi, Theta\) = \(0.0, 61.0\) of key \(0.0, 135.0\) is not a node of the declared grid "
     r"Phi 0.0 to 30.0 in 2 nodes x Theta 60.0 to 120.0 in 3 nodes"),
    (dict(exit_yz=((-0.001, -5e-6), (-0.001, 5e-6), (0.0, -5e-6), (0.0, -5e-6), (0.001, -5e-6), (0.001, 5e-6))),
     r"exit point \(Y, Z\) = \(0.0, -5e-06\) m appears twice in key \(0.0, 135.0\)"),
    (dict(exit_yz=((-0.001, 0.0), (0.0, 0.0), (0.0015, 0.0))),
     r"exit point \(Y, Z\) = \(0.0, 0.0\) m of key \(0.0, 135.0\) is not a node of the declared lattice"),
], ids=["duplicate-far-field-point", "uneven-far-field-axis", "duplicate-exit-point", "uneven-exit-axis"])
def test_c3_c4_per_key_grids(tmp_path, kw, match):
    """CDX-022: the external review's fixtures (and an uneven exit axis) pass the dataset-wide counts; the per-key
    node checks refuse them."""
    base = str(tmp_path)
    stem = write_grid(base, **kw)
    sidecar.write(base, stem, tiny_sidecar(base, stem))
    with pytest.raises(ValueError, match=f"^BBR012: c_500GHz_Ephi=0: {match}"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_c3_keys_differing_in_compensating_ways(tmp_path):
    """CDX-022: (0, 135) holds (0, 60) and (30, 120), (0, 180) holds (30, 60) and (0, 120): the union
    is the declared 2 x 2 grid, so only the per-key check sees that each key holds half of it."""
    base = str(tmp_path)
    stem = write_grid(base, ff=((0, 60), (30, 120)))
    for e in (0, 1):
        rewrite(base, stem, e, "far_field.csv",
                lambda t: t.replace(",0.0,180.0,0,60,", ",0.0,180.0,30,60,").replace(",0.0,180.0,30,120,", ",0.0,180.0,0,120,"))
    sidecar.write(base, stem, tiny_sidecar(base, stem))
    with pytest.raises(ValueError, match=r"^BBR012: c_500GHz_Ephi=0: key \(0.0, 135.0\) holds 2 of the 4 nodes"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_check_full_refuses_a_transmitting_key_without_support(tmp_path):
    """CDX-008 in the preflight: Ephi=0 transmits (T0 = 0.25 at (0, 135)) but its far field is zero
    everywhere; BBRHFSSData refuses it at load (BBR027), so check_full must too, before a run."""
    base = str(tmp_path)
    stem = write_grid(base)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 0, "far_field.csv", lambda t: t.replace(",0.1,0.0,0.2,0.05\n", ",0.0,0.0,0.0,0.0\n"))
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR027: c_500GHz key \(0.0, 135.0\): T0 = 0.25 > 0, but its far field "
                                         r"has no outward intensity"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_c3_pole_rows_are_valid(tmp_path):
    """Theta = 0 at two Phi is one direction but two grid nodes: valid (CDX-022 compares coordinate
    pairs, not directions)."""
    base = str(tmp_path)
    stem = write_grid(base, ff=((0, 0), (90, 0), (0, 90), (90, 90)))
    sidecar.write(base, stem, tiny_sidecar(base, stem))
    sidecar.check_full(base, "c", stem, 500.0)


def test_c4_omitted_requires_every_inside_node(tmp_path):
    """With "omitted" the per-key exit set is the lattice inside the section: R = 15 um at a 10 um step
    holds -10, 0 and +10 um; dropping the centre point is refused."""
    base = str(tmp_path)
    stem = write_disc(base)
    for e in (0, 1):
        rewrite(base, stem, e, "waveguide.csv", lambda t: "".join(l for l in t.splitlines(keepends=True)
                                                                if ",0,0.0," not in l))
    sc = disc_sidecar(base, stem, "omitted")
    sc["exit_field"]["grid"]["y_e"] = {"min": -2e-5, "max": 2e-5, "count": 5}
    with pytest.raises(ValueError, match=r"^BBR012: c_500GHz_Ephi=0: key \(0.0, 180.0\) holds 2 exit points; every "
                                         r"key holds each of the 3 declared lattice nodes inside the section once"):
        sidecar.check_csvs(base, stem, sc)


def swap_pairs(text):
    """text with its data rows traded in consecutive pairs (each key's two rows swap places)."""
    head, *rows = text.splitlines(keepends=True)
    for i in range(0, len(rows) - 1, 2):
        rows[i], rows[i + 1] = rows[i + 1], rows[i]
    return head + "".join(rows)


@pytest.mark.parametrize("name,kw,match", [
    ("waveguide.csv", dict(exit_yz=((-0.001, 0.0), (0.001, 0.0))),
     r"exit point 0 lies at \(X, Y, Z\) = \(0.0, 0.001, 0.0\), the Ephi=0 one at \(0.0, -0.001, 0.0\)"),
    ("far_field.csv", dict(),
     r"far-field row 0 lies at \(Phi, Theta\) = \(0.0, 120.0\), the Ephi=0 one at \(0.0, 60.0\)"),
], ids=["exit-rows-swapped", "far-field-rows-swapped"])
def test_check_csvs_pairs_rows_by_order(tmp_path, name, kw, match):
    """CDX-002: Ephi=1 with each key's two rows swapped keeps every key, grid, count and checksum, but
    BBRHFSSData pairs the rows by order (BBR012); the preflight must refuse it too."""
    base = str(tmp_path)
    stem = write_grid(base, **kw)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 1, name, swap_pairs)
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=rf"^BBR012: c_500GHz: key \(0.0, 135.0\): Ephi=1 {match}"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_check_csvs_unequal_row_counts(tmp_path):
    """CDX-002's unequal-count case: one Ephi=1 key with a row fewer."""
    base = str(tmp_path)
    stem = write_grid(base, exit_yz=((-0.001, 0.0), (0.001, 0.0)))
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 1, "waveguide.csv", lambda t: t.replace(t.splitlines(keepends=True)[1], "", 1))
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR012: c_500GHz_Ephi=1: exit-point rows per incidence key differ"):
        sidecar.check_full(base, "c", stem, 500.0)


@pytest.mark.parametrize("name,old,new,match", [
    ("far_field.csv", ",0.2,0.05\n", ",inf,0.05\n", r"line 2, column rEtheta_real: 'inf' is not a finite number"),
    ("far_field.csv", ",0.2,0.05\n", ",-inf,0.05\n", r"line 2, column rEtheta_real: '-inf' is not a finite number"),
    ("far_field.csv", ",0.2,0.05\n", ",nan,0.05\n", r"empty field or NaN \(line 2, column rEtheta_real\)"),
    ("waveguide.csv", ",0.25,1.0,", ",inf,1.0,", r"line 2, column OutgoingPower: 'inf' is not a finite number"),
    ("waveguide.csv", ",0.25,1.0,", ",abc,1.0,", r"line 2, column OutgoingPower: 'abc' is not a finite number"),
], ids=["inf", "minus-inf", "nan", "outgoing-power-inf", "not-a-number"])
def test_check_full_rejects_non_finite(tmp_path, name, old, new, match):
    """CDX-001: a value BBRHFSSData would refuse at load (BBR013) fails the preflight too, with the
    file, line and column, even with fresh checksums."""
    base = str(tmp_path)
    stem = write_grid(base)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 1, name, lambda t: t.replace(old, new, 1))
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=rf"^BBR013: .*c_500GHz_Ephi=1/{name}: .*{match}"):
        sidecar.check_full(base, "c", stem, 500.0)


def test_check_full_rejects_a_foreign_freq_token(tmp_path):
    """BBRHFSSData refuses a Freq token more than 0.1 % from the directory frequency (BBR009)."""
    base = str(tmp_path)
    stem = write_grid(base)
    sc = tiny_sidecar(base, stem)
    rewrite(base, stem, 1, "far_field.csv", lambda t: t.replace("500GHz,1,0.0,180.0,", "450GHz,1,0.0,180.0,", 1))
    write_with_files(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR009: c_500GHz_Ephi=1/far_field.csv: Freq column '450GHz'"):
        sidecar.check_full(base, "c", stem, 500.0)


# --- fix round 2: the physics-only invariant, the disc modes, load(), malformed CSVs -------------

def built(tmp_path):
    base = str(tmp_path)
    return tiny_sidecar(base, write_tiny(base))


def test_invariant_keys(tmp_path):
    assert list(sidecar.invariant(built(tmp_path))) == ["schema", "excitation", "frames", "symmetry", "boundaries",
                                                       "geometry", "modes", "cross_section"]


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


LEGACY = ["InfParallelPlate_crack1Rohan_500GHz", "InfParallelPlate_crack2_500GHz"]


@pytest.mark.parametrize("stem", LEGACY)
def test_legacy_sidecar_passes_full_check(data_root, stem):
    base = os.path.join(data_root, "waveguides")
    sc = sidecar.check_full(base, stem.rsplit("_", 1)[0], stem, 500.0)
    assert sc["modes"]["mode"] == "TE10" and sc["modes"]["propagating_count"] == 33
    assert sc["frames"]["pose_rule"] == "legacy" and sc["provenance"]["inferred"] == ["boundaries.walls"]


# The HFSS writer's frequency-independent values for the two cracks (Blackbody-Simulations run_job on the
# crack base geometries, 2026-10-05). A frequency from that writer must agree with the legacy 500 GHz
# sidecar (same_invariant), or a tree that mixes them is refused (BBR024).
HFSS_WRITER = {
    "InfParallelPlate_crack1Rohan": dict(extent={"g": 0.05, "l": 10.0, "p": 1.0}, z_half=2.5e-05,
                                         filter=2997.92458, onsets=7, count=13970),
    "InfParallelPlate_crack2": dict(extent={"g": 0.1, "l": 10.0, "p": 1.5}, z_half=5e-05,
                                    filter=1498.96229, onsets=14, count=28061),
}


@pytest.mark.parametrize("stem", LEGACY)
def test_legacy_invariant_matches_hfss_writer(data_root, stem):
    w = HFSS_WRITER[stem.rsplit("_", 1)[0]]
    legacy = sidecar.load(os.path.join(data_root, "waveguides"), stem)
    onsets = [{"n": 0, "mode": "TE10", "cutoff_ghz": 14.9896229}] + [
        {"n": k, "mode": sidecar.mode_label("TE", 0, k), "cutoff_ghz": k * w["filter"]} for k in range(1, w["onsets"])]
    writer = {
        "frames": {"hfss_global_axes_in_canonical": {"x": [0.0, 0.0, -1.0], "y": [0.0, 1.0, 0.0], "z": [1.0, 0.0, 0.0]},
                   "exit_cs_axes_in_canonical": {"x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0], "z": [0.0, 0.0, 1.0]}},
        "symmetry": {"end_to_end": True, "mirror_g": True, "mirror_l": True, "rotational": False},
        "boundaries": {"entrance": "radiation", "exit": "radiation", "walls": "PEC"},
        "geometry": {"shape": "box", "extent_mm": w["extent"]},
        "exit_field": {"cross_section": {"shape": "rectangle", "y_e_half_m": 0.005, "z_e_half_m": w["z_half"]}},
        "modes": {"mode": "TE10", "cutoff_ghz": 14.9896229, "list_limit_ghz": 20000.0,
                  "polarization_filter_limit_ghz": w["filter"], "gap_family_onsets": onsets,
                  "mode_count_below_limit": w["count"], "propagating_count": 33, "basis": "the HFSS writer's wording"},
    }
    assert sidecar.same_invariant(legacy, writer)
    assert legacy["frames"]["canonical"] == "p,l,g; p x l = g; Geant4 crack-local (x,y,z) = (p,l,g)"


def test_legacy_sidecars_validator(repo_root, data_root):
    import subprocess, sys
    r = subprocess.run([sys.executable, str(repo_root / "validation/check_dataset_sidecars.py"),
                        os.path.join(data_root, "waveguides")], capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip().endswith("RESULT: PASS"), r.stdout + r.stderr
    assert "2 dataset(s) checked, 0 failing" in r.stdout.splitlines(), r.stdout
    base = os.path.join(data_root, "waveguides")
    assert [l for l in r.stdout.splitlines() if l.startswith("  ")] == [f"  PASS {base}/{s}" for s in LEGACY]


# --- invariant_diff and the validator on tiny trees ------------------------------

@pytest.mark.parametrize("edit,blocks", [
    (lambda s: s["geometry"]["extent_mm"].update(p=4.0), ["geometry"]),
    (lambda s: s["symmetry"].update(rotational=True), ["symmetry"]),
    (lambda s: s["boundaries"].update(walls="finite-conductivity"), ["boundaries"]),
    (lambda s: s["modes"].update(cutoff_ghz=s["modes"]["cutoff_ghz"] * (1 + 1e-6)), ["modes"]),
    (lambda s: (s["boundaries"].update(walls="x"), s["exit_field"]["cross_section"].update(y_e_half_m=2.5e-3)),
     ["boundaries", "cross_section"]),
    (lambda s: s["frames"].update(pose_rule="legacy"), []),
    (lambda s: s["modes"].update(propagating_count=1), []),
], ids=["geometry", "symmetry", "boundaries", "modes", "two-blocks", "pose-rule", "propagating-count"])
def test_invariant_diff_names_blocks(tmp_path, edit, blocks):
    a = built(tmp_path)
    b = copy.deepcopy(a)
    edit(b)
    assert sidecar.invariant_diff(a, b) == blocks and sidecar.invariant_diff(b, a) == blocks
    assert sidecar.same_invariant(a, b) == (not blocks)


def tiny_dataset(base, label="500GHz", dataset_id="c", edit=None):
    """A complete tiny dataset (both Ephi directories and a passing sidecar) under base."""
    os.makedirs(base, exist_ok=True)
    stem = write_tiny(base, dataset_id=dataset_id, label=label)
    sc = tiny_sidecar(base, stem, dataset_id=dataset_id, label=label)
    if edit:
        edit(sc)
    sidecar.write(base, stem, sc)
    return stem


def run_validator(repo_root, *dirs):
    import subprocess, sys
    r = subprocess.run([sys.executable, str(repo_root / "validation/check_dataset_sidecars.py"), *map(str, dirs)],
                       capture_output=True, text=True)
    assert not r.stderr, r.stderr                  # no traceback, whatever the input
    lines = r.stdout.splitlines()
    return r.returncode, [l for l in lines if l.startswith("  ")], lines[-2:]


def test_validator_tiny_tree_passes_and_ignores_other_files(tmp_path, repo_root):
    base = tmp_path / "wg"
    tiny_dataset(str(base))
    tiny_dataset(str(base), label="600GHz")
    (base / "SHA256SUMS").write_text("x\n")
    (base / "manifest.json").write_text("{}\n")
    (base / "notes_Ephi=0").write_text("a file, not a directory\n")
    rc, rows, tail = run_validator(repo_root, base)
    assert rc == 0 and rows == [f"  PASS {base}/c_500GHz", f"  PASS {base}/c_600GHz"]
    assert tail == ["2 dataset(s) checked, 0 failing", "RESULT: PASS"]


@pytest.mark.parametrize("orphan,missing", [
    (lambda b: os.makedirs(b / "c_600GHz_Ephi=1"), "c_600GHz_Ephi=0/, c_600GHz.dataset.json"),
    (lambda b: os.makedirs(b / "c_600GHz_Ephi=0"), "c_600GHz_Ephi=1/, c_600GHz.dataset.json"),
    (lambda b: (b / "c_600GHz.dataset.json").write_text("{}\n"), "c_600GHz_Ephi=0/, c_600GHz_Ephi=1/"),
], ids=["orphan-ephi1", "orphan-ephi0", "orphan-sidecar"])
def test_validator_fails_an_incomplete_dataset(tmp_path, repo_root, orphan, missing):
    base = tmp_path / "wg"
    tiny_dataset(str(base))
    orphan(base)
    rc, rows, tail = run_validator(repo_root, base)
    assert rc == 1 and rows == [f"  PASS {base}/c_500GHz",
                                f"  FAIL {base}/c_600GHz: ValueError: incomplete dataset: missing {missing}"]
    assert tail == ["2 dataset(s) checked, 1 failing", "RESULT: FAIL"]


def test_validator_fails_a_missing_sidecar(tmp_path, repo_root):
    base = tmp_path / "wg"
    tiny_dataset(str(base))
    write_tiny(str(base), label="600GHz")
    rc, rows, _ = run_validator(repo_root, base)
    assert rc == 1 and rows[1] == f"  FAIL {base}/c_600GHz: ValueError: incomplete dataset: missing c_600GHz.dataset.json"


def test_validator_cross_frequency_disagreement_names_the_block(tmp_path, repo_root):
    base = tmp_path / "wg"
    tiny_dataset(str(base))
    tiny_dataset(str(base), label="600GHz", edit=lambda s: s["symmetry"].update(rotational=True))
    rc, rows, tail = run_validator(repo_root, base)
    assert rc == 1 and rows == [
        f"  PASS {base}/c_500GHz",
        f"  FAIL {base}/c_600GHz: disagrees with c_500GHz on the frequency-independent block(s) symmetry"]
    assert tail == ["2 dataset(s) checked, 1 failing", "RESULT: FAIL"]


def test_validator_empty_or_missing_directory_argument(tmp_path, repo_root):
    good, empty, absent = tmp_path / "wg", tmp_path / "empty", tmp_path / "absent"
    tiny_dataset(str(good))
    empty.mkdir()
    rc, rows, tail = run_validator(repo_root, good, empty, absent)
    assert rc == 1 and rows[:2] == [f"  PASS {good}/c_500GHz", f"  FAIL {empty}: no dataset found"]
    assert rows[2].startswith(f"  FAIL {absent}: cannot list the directory (FileNotFoundError: ") and len(rows) == 3
    assert tail == ["1 dataset(s) checked, 0 failing; 2 directory argument(s) failing", "RESULT: FAIL"]


def test_validator_malformed_stem_and_bad_sidecar_do_not_hide_the_rest(tmp_path, repo_root):
    base = tmp_path / "wg"
    for part in ("a_Ephi=0", "a_Ephi=1"):
        (base / part).mkdir(parents=True)
    (base / "a.dataset.json").write_text("{}\n")
    stem = tiny_dataset(str(base), dataset_id="b")
    (base / f"{stem}.dataset.json").write_text("[1, 2\n")        # not JSON
    tiny_dataset(str(base))
    rc, rows, tail = run_validator(repo_root, base)
    assert rc == 1 and rows[0] == f"  FAIL {base}/a: ValueError: stem 'a' is not <id>_<freq>GHz"
    assert rows[1].startswith(f"  FAIL {base}/b_500GHz: ValueError: BBR024: ") and "not valid JSON" in rows[1]
    assert rows[2] == f"  PASS {base}/c_500GHz" and len(rows) == 3
    assert tail == ["3 dataset(s) checked, 2 failing", "RESULT: FAIL"]


def test_validator_fails_a_tampered_checksum(tmp_path, repo_root):
    base = tmp_path / "wg"
    stem = tiny_dataset(str(base))
    with open(base / f"{stem}_Ephi=0" / "waveguide.csv", "a") as fh:
        fh.write("\n")
    rc, rows, tail = run_validator(repo_root, base)
    assert rc == 1 and rows[0].startswith(f"  FAIL {base}/c_500GHz: ValueError: BBR024: ") and "sha256" in rows[0]
    assert tail[-1] == "RESULT: FAIL"


# --- schema 2.x (Tasks 1-3 of the reader plan: OR-028, OR-031, OR-034, OR-036, OR-038, HF-040, HF-047, HF-049) ----
# The reference is tests/schema2_fixture.py, written independently of bbrsim.sidecar and the same content
# as tests/HFSSFixture.hh (NineKeyGrid2, SidecarJson2).

import hashlib

import schema2_fixture as fx

ID = "crack"


def setp(d, path, value):
    """d[a][b]...[z] = value for a slash path; a numeric step indexes a list."""
    *head, last = path.split("/")
    for k in head:
        d = d[int(k)] if isinstance(d, list) else d[k]
    if isinstance(d, list):
        d[int(last)] = value
    else:
        d[last] = value


def mut2(fn, **kw):
    s = fx.sidecar2(**kw)
    fn(s)
    return s


def check2(sc, stem=fx.STEM, freq=500.0):
    sidecar.check(sc, ID, stem, freq)


def s2_dataset(tmp_path, **kw):
    base = str(tmp_path)
    return base, fx.write_nine_key_grid2(base, **kw)


def edit_key(path, phi, theta, fn):
    """Apply fn(fields, n) to every data row of key (phi, theta) in a CSV (n counts the key's rows)."""
    lines = open(path).read().splitlines()
    n = 0
    for i, line in enumerate(lines[1:], start=1):
        v = line.split(",")
        if float(v[2]) == phi and float(v[3]) == theta:
            fn(v, n)
            n += 1
            lines[i] = ",".join(v)
    assert n, "no row of the key"
    open(path, "w").write("\n".join(lines) + "\n")


def wg_path(base, stem, e):
    return os.path.join(base, f"{stem}_Ephi={e}", "waveguide.csv")


def test_s2_parses(tmp_path):
    sc = fx.sidecar2()
    check2(sc)
    assert sidecar.schema_major(sc) == 2 and sidecar.schema_major(fx.sidecar2(stem="crack_500GHz")) == 2
    recs = sidecar.key_records(sc)
    assert len(recs) == 9
    assert recs[0]["phi"] == 0 and recs[0]["theta"] == 90 and recs[1]["theta"] == 135 and recs[3]["phi"] == 45
    dark = [r for r in recs if r["phi"] == 0 and r["theta"] == 180][0]
    assert dark["dark1"] and not dark["dark0"] and dark["t_flux0"] == 0.6 and dark["t_flux1"] == 0.0
    assert dark["c_k"] == pytest.approx(1.0, abs=1e-12) and dark["w1"] == 1.0
    a = [r for r in recs if (r["phi"], r["theta"]) == (45, 135)][0]
    assert not a["dark0"] and not a["dark1"] and a["w0"] == 0.9 and a["t_flux1"] == 0.3
    assert a["c_k"] == pytest.approx(math.cos(math.radians(45)), rel=1e-12)
    e = [r for r in recs if (r["phi"], r["theta"]) == (90, 90)][0]
    assert 0.083 < e["c_k"] < 0.085
    # the disc variant (R2: uniform_disc, TE11 lowest) parses too, and so does the 1.x schema untouched
    check2(fx.sidecar2(disc=True))
    assert sidecar.schema_major({"schema_version": "1.0"}) == 1


def test_s2_version_gate():
    for v in ("3.0", "0.9", "20.0", "12", "2x"):
        with pytest.raises(ValueError, match=r"^BBR024: .*BBRsim reads 1\.x and 2\.x"):
            check2(mut2(lambda s: s.update(schema_version=v)))
    with pytest.raises(ValueError, match=r"^BBR024"):
        sidecar.schema_major({"schema_version": 2})
    with pytest.raises(ValueError, match=r"^BBR024"):
        sidecar.schema_major({})
    check2(mut2(lambda s: s.update(schema_version="2.7", x_future=[1, 2])))
    check2(mut2(lambda s: s.update(schema_version="2")))
    assert sidecar.schema_major({"schema_version": "1"}) == 1 and sidecar.schema_major({"schema_version": "1.4"}) == 1


S2_STRING_EDITS = [
    ("model/kind", "slot_in_wall", "model.kind"),
    ("transmittance/definition", "OutgoingPower / IngoingPower", "transmittance.definition"),
    ("transmittance/numerator", "exit_abs_re_s", "transmittance.numerator"),
    ("transmittance/incoming", "area_times_intensity_times_cos", "transmittance.incoming"),
    ("transmittance/above_one", "load_time_rescale", "transmittance.above_one"),
    ("exit_field/tangential_components", ["Ex", "Ez"], "exit_field.tangential_components"),
    ("exit_field/quadrature", "uniform_disc", "exit_field.quadrature"),
    ("exit_field/sampler", "abs_e_squared", "exit_field.sampler"),
    ("exit_field/scale/rule", "none", "exit_field.scale.rule"),
    ("cross_term/density", "all_components", "cross_term.density"),
    ("cross_term/above_limit", "ignore", "cross_term.above_limit"),
    ("far_field/normalization", "absolute", "far_field.normalization"),
    ("exit_field/column_order_check/rule", "none", "exit_field.column_order_check.rule"),
    ("exit_field/column_order_check/reference/passed", False, "exit_field.column_order_check.reference.passed"),
    ("exit_field/plane_tolerance_m", 1e-9, "exit_field.plane_tolerance_m"),
    ("exit_field/columns", ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower", "X", "Y", "Z",
                            "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag", "InterceptFactor"],
     "exit_field.columns"),
]


def test_s2_fixed_strings():
    # the module's dict carries every fixed string of the field list
    assert sidecar.SCHEMA2_STRINGS["transmittance.definition"] == "OutgoingPower / (IngoingPower * InterceptFactor)"
    assert sidecar.SCHEMA2_STRINGS["exit_field.scale.rule"] == \
        "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w"
    assert sidecar.WG_COLUMNS_2[9:15] == ["Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag"]
    assert sidecar.WG_COLUMNS_2[15] == "InterceptFactor" and len(sidecar.WG_COLUMNS_2) == 16
    for path, value, field in S2_STRING_EDITS:
        with pytest.raises(ValueError, match=r"^BBR025: .*" + field.replace(".", r"\.")):
            check2(mut2(lambda s: setp(s, path, value)))
    # trapezoid is the only quadrature of a rectangle, and uniform_disc of a disc
    with pytest.raises(ValueError, match=r"^BBR025: .*exit_field\.quadrature"):
        check2(mut2(lambda s: s["exit_field"].update(quadrature="trapezoid"), disc=True))
    # a missing block is BBR024, a mistyped string too
    with pytest.raises(ValueError, match=r"^BBR024: missing field model"):
        check2(mut2(lambda s: s.pop("model")))
    with pytest.raises(ValueError, match=r"^BBR024: transmittance\.numerator"):
        check2(mut2(lambda s: setp(s, "transmittance/numerator", 3)))
    with pytest.raises(ValueError, match=r"^BBR024: missing field exit_field\.scale"):
        check2(mut2(lambda s: s["exit_field"].pop("scale")))
    with pytest.raises(ValueError, match=r"^BBR024: conventions"):
        check2(mut2(lambda s: s.update(conventions=[1])))
    # the 1.x column list is not the 2.x one, and a 2.x sidecar naming the 1.0 text is refused
    with pytest.raises(ValueError, match=r"^BBR025: .*exit_field\.columns"):
        check2(mut2(lambda s: s["exit_field"].update(columns=list(sidecar.WG_COLUMNS))))


def test_s2_exit_cs_field_frame_refused():
    """CDX-035: 2.x takes E_t as the HFSS-global (Ex, Ey), so exit-CS field components are refused even when
    declared consistently."""
    with pytest.raises(ValueError, match=r"^BBR025: exit_field\.field_components_frame is exit_cs"):
        check2(mut2(lambda s: s["exit_field"].update(field_in_ref_cs=True, field_components_frame="exit_cs")))
    with pytest.raises(ValueError, match=r"^BBR025: exit_field\.field_in_ref_cs and field_components_frame contradict"):
        check2(mut2(lambda s: s["exit_field"].update(field_components_frame="exit_cs")))


def test_s1_exit_cs_field_frame_accepted(tmp_path):
    """Under 1.x exit-CS components declared consistently stay accepted: |E|^2 and the overlap do not depend
    on the frame (CDX-035 leaves schema 1 as it was)."""
    base = str(tmp_path)
    stem = write_tiny(base)
    sc = tiny_sidecar(base, stem)
    sc["exit_field"].update(field_in_ref_cs=True, field_components_frame="exit_cs")
    sidecar.check(sc, "c", stem, 500.0)


S2_LISTS = [("transmittance/c_key", "c_k", False), ("transmittance/per_key", "polarizations/0/t_flux", True),
            ("exit_field/scale/per_key", "w0", False), ("far_field/ff_over_flux", "polarizations/0", True)]


@pytest.mark.parametrize("lst,number,pols", S2_LISTS, ids=[x[0] for x in S2_LISTS])
def test_s2_per_key_lists(lst, number, pols):
    def get(s):
        d = s
        for k in lst.split("/"):
            d = d[k]
        return d

    # a duplicate key, a missing key, a key off the declared grid by 1e-6, a string angle, an object
    with pytest.raises(ValueError, match=r"^BBR024: .*duplicate"):
        check2(mut2(lambda s: get(s).append(copy.deepcopy(get(s)[0]))))
    with pytest.raises(ValueError, match=r"^BBR024: .*no record for declared key"):
        check2(mut2(lambda s: get(s).pop()))
    with pytest.raises(ValueError, match=r"^BBR024: .*does not match a declared key"):
        check2(mut2(lambda s: get(s)[0].update(theta_deg=90.0 + 1e-6)))
    with pytest.raises(ValueError, match=r"^BBR024: .*phi_deg"):
        check2(mut2(lambda s: get(s)[0].update(phi_deg="0")))
    with pytest.raises(ValueError, match=r"^BBR024: .*must be a list"):
        check2(mut2(lambda s: setp(s, lst, {})))
    if pols:
        with pytest.raises(ValueError, match=r"^BBR024: .*polarizations"):
            check2(mut2(lambda s: get(s)[0]["polarizations"].pop("1")))
        with pytest.raises(ValueError, match=r"^BBR024: .*polarizations"):
            check2(mut2(lambda s: get(s)[0].pop("polarizations")))
        with pytest.raises(ValueError, match=r"^BBR024: .*polarizations"):
            check2(mut2(lambda s: get(s)[0]["polarizations"].update({"2": 1.0})))
    # a non-finite or non-numeric number is refused (NaN and an overflowed 1e999 read back as non-finite)
    for bad in (float("nan"), float("inf"), "nan"):
        with pytest.raises(ValueError, match=r"^BBR024"):
            check2(mut2(lambda s: setp(s, f"{lst}/0/{number}", bad)))
    # 5e-10 off the declared key still matches (1e-9)
    s = fx.sidecar2()
    get(s)[0]["theta_deg"] = 90.0 + 5e-10
    get(s)[3]["phi_deg"] = 45.0 - 5e-10
    check2(s)


def test_s2_per_key_other_checks():
    # ff_over_flux: null on a lit basis is refused, null on a dark basis parses (the fixture has one)
    with pytest.raises(ValueError, match=r"^BBR024: .*null"):
        check2(mut2(lambda s: setp(s, "far_field/ff_over_flux/0/polarizations/0", None)))
    assert fx.sidecar2()["far_field"]["ff_over_flux"][2]["polarizations"]["1"] is None
    check2(mut2(lambda s: setp(s, "far_field/ff_over_flux/0/polarizations/0", 0.0)))
    with pytest.raises(ValueError, match=r"^BBR024: .*finite number"):
        check2(mut2(lambda s: setp(s, "far_field/ff_over_flux/0/polarizations/0", "x")))
    # a scale weight must be positive
    with pytest.raises(ValueError, match=r"^BBR024: .*w0"):
        check2(mut2(lambda s: setp(s, "exit_field/scale/per_key/0/w0", 0.0)))
    # per_key c_k disagreeing with c_key is a self-contradiction: BBR024
    with pytest.raises(ValueError, match=r"^BBR024: .*c_k"):
        check2(mut2(lambda s: setp(s, "transmittance/per_key/1/c_k", 0.5)))


def test_s2_cell90_values():
    assert sidecar.cell_90([0, 45, 90, 135, 180]) == pytest.approx(0.08384881979, abs=1e-10)
    assert sidecar.cell_90([90, 105, 120, 135, 150, 165, 180]) == pytest.approx(0.03663305469, abs=1e-10)
    assert sidecar.cell_90([0, 45, 90, 135, 180]) == pytest.approx(0.0838488197937161, abs=1e-14)
    assert sidecar.cell_90([90, 105, 120, 135, 150, 165, 180]) == pytest.approx(0.036633054690581604, abs=1e-14)
    assert sidecar.cell_90([120, 180]) is None
    assert sidecar.cell_90([0, 45, 90]) is None                       # nothing above 90 deg
    assert sidecar.intercept_factor(135.0, [90, 135, 180]) == pytest.approx(math.cos(math.radians(45)), rel=1e-14)
    assert sidecar.intercept_factor(180.0, [90, 135, 180]) == pytest.approx(1.0, abs=1e-15)
    assert sidecar.intercept_factor(90.0, [90, 135, 180]) == sidecar.cell_90([90, 135, 180])
    assert sidecar.intercept_factor(45.0, [0, 45, 90, 135, 180]) == pytest.approx(math.cos(math.radians(45)), rel=1e-14)
    # c_key at (0, 135) off by 2e-6 relative is refused (BBR025), off by 5e-7 parses
    c135 = math.cos(math.radians(45))
    with pytest.raises(ValueError, match=r"^BBR025: .*c_k"):
        check2(mut2(lambda s: setp(s, "transmittance/c_key/1/c_k", c135 * (1 + 2e-6))))
    check2(mut2(lambda s: setp(s, "transmittance/c_key/1/c_k", c135 * (1 + 5e-7))))
    # cell_90 and the c_k at 90 deg must follow the declared grid
    with pytest.raises(ValueError, match=r"^BBR025: .*cell_90"):
        check2(mut2(lambda s: setp(s, "transmittance/cell_90", 0.08384881979 * 1.001)))
    with pytest.raises(ValueError, match=r"^BBR025: .*c_k"):
        check2(mut2(lambda s: setp(s, "transmittance/c_key/0/c_k", 0.19134)))
    check2(mut2(lambda s: setp(s, "transmittance/cell_90", 0.0838488197937161 * (1 + 5e-7))))
    # the 35-key style grid: Delta = 15 deg gives the 0.0366 cell
    th = (90, 105, 120, 135, 150, 165, 180)
    s = fx.sidecar2(thetas=th)
    check2(s)
    assert s["transmittance"]["cell_90"] == pytest.approx(0.03663305469, abs=1e-10) and len(sidecar.key_records(s)) == 21


def test_s2_no_theta_90_refused():
    # a declared grid without 90 deg cannot give the cell: BBR025, before any per-key matching
    s = fx.sidecar2()
    s["excitation"]["incident_theta_deg"] = [135.0, 180.0]
    with pytest.raises(ValueError, match=r"^BBR025: .*theta = 90"):
        check2(s)
    s["excitation"]["incident_theta_deg"] = [90.0]                    # nothing above 90 deg either
    with pytest.raises(ValueError, match=r"^BBR025: .*theta = 90"):
        check2(s)


def test_s2_dark_tflux_bbr024():
    # OR-038: a dark basis has t_flux exactly 0; the sidecar contradicting itself is BBR024
    bad = mut2(lambda s: setp(s, "transmittance/per_key/2/polarizations/1/t_flux", 1e-15))
    with pytest.raises(ValueError, match=r"^BBR024: .*key \(phi 0, theta 180\).*polarization 1.*1e-15.*OR-038"):
        check2(bad)
    with pytest.raises(ValueError, match=r"^BBR024: .*polarization 0"):
        check2(mut2(lambda s: setp(s, "exit_field/scale/per_key/0/dark0", True)))
    check2(fx.sidecar2())                                              # the fixture's own dark key parses
    # a lit basis with t_flux 0 is fine
    check2(mut2(lambda s: setp(s, "transmittance/per_key/0/polarizations/0/t_flux", 0.0)))


def test_s2_emitter():
    base = "transmittance/cell_90_assumes"
    for path, value, field in ((f"{base}/emitter", "lambertian", "emitter"), (f"{base}/cos_floor", 0.02, "cos_floor"),
                               (f"{base}/band_deg", [60.0, 90.0], "band_deg")):
        with pytest.raises(ValueError, match=r"^BBR025: .*cell_90_assumes\." + field):
            check2(mut2(lambda s: setp(s, path, value)))
    with pytest.raises(ValueError, match=r"^BBR024: .*band_deg"):
        check2(mut2(lambda s: s["transmittance"]["cell_90_assumes"].pop("band_deg")))
    with pytest.raises(ValueError, match=r"^BBR024: .*band_deg"):
        check2(mut2(lambda s: setp(s, f"{base}/band_deg", [67.5])))
    check2(mut2(lambda s: setp(s, f"{base}/band_deg", [67.5 + 5e-10, 90.0])))


def test_s2_grid_step():
    for ax, bad in (("y_e", 1.001e-3), ("z_e", 2.4e-5)):
        with pytest.raises(ValueError, match=rf"^BBR024: .*exit_field\.grid\.{ax}\.step"):
            check2(mut2(lambda s: setp(s, f"exit_field/grid/{ax}/step", bad)))
        with pytest.raises(ValueError, match=rf"^BBR024: .*exit_field\.grid\.{ax}\.step"):
            check2(mut2(lambda s: s["exit_field"]["grid"][ax].pop("step")))
        # choice 3: a rectangle axis with one node has no step and no trapezoid
        with pytest.raises(ValueError, match=rf"^BBR024: .*exit_field\.grid\.{ax} needs count >= 2"):
            check2(mut2(lambda s: setp(s, f"exit_field/grid/{ax}/count", 1)))
    with pytest.raises(ValueError, match=r"^BBR024: .*step must be > 0"):
        check2(mut2(lambda s: setp(s, "exit_field/grid/y_e/step", 0.0)))
    check2(mut2(lambda s: setp(s, "exit_field/grid/y_e/step", 1e-3 * (1 + 5e-7))))


def test_s2_regime_cutoffs():
    r1 = sidecar.regime_cutoffs({"shape": "rectangle", "y_e_half_m": 5e-3, "z_e_half_m": 25e-6})
    assert len(r1) == 7 and [r1[0]["mode"], r1[1]["mode"], r1[6]["mode"]] == ["TE10", "TE01", "TE06"]
    assert r1[0]["cutoff_ghz"] == pytest.approx(14.9896229, rel=1e-6)
    assert r1[1]["cutoff_ghz"] == pytest.approx(2997.92458, rel=1e-6)
    assert r1[6]["cutoff_ghz"] == pytest.approx(17987.5475, rel=1e-6)
    r2 = sidecar.regime_cutoffs({"shape": "rectangle", "y_e_half_m": 5e-3, "z_e_half_m": 50e-6})
    assert len(r2) == 14 and (r2[0]["mode"], r2[9]["mode"], r2[10]["mode"], r2[13]["mode"]) == \
        ("TE10", "TE09", "TE0,10", "TE0,13")
    assert r2[13]["cutoff_ghz"] == pytest.approx(19486.5098, rel=1e-6)
    rd = sidecar.regime_cutoffs({"shape": "disc", "radius_m": 50e-6})
    assert [e["mode"] for e in rd] == ["TE11", "TM01", "TE21"]
    assert [e["cutoff_ghz"] for e in rd] == pytest.approx([1756.98, 2294.85, 2914.56], rel=1e-5)
    # the fixture's own list is accepted; an entry off by 2e-6 is refused (BBR025), by 5e-7 it parses
    s = fx.sidecar2()
    assert len(s["modes"]["regime_cutoffs_ghz"]) == 7
    c1 = s["modes"]["regime_cutoffs_ghz"][1]["cutoff_ghz"]
    for tag in ("regime_cutoffs_ghz", "TE01"):
        with pytest.raises(ValueError, match=r"^BBR025: .*" + tag):
            check2(mut2(lambda s: setp(s, "modes/regime_cutoffs_ghz/1/cutoff_ghz", c1 * (1 + 2e-6))))
    check2(mut2(lambda s: setp(s, "modes/regime_cutoffs_ghz/1/cutoff_ghz", c1 * (1 + 5e-7))))
    # a missing, an extra and a relabelled entry are refused
    with pytest.raises(ValueError, match=r"^BBR025: .*regime_cutoffs_ghz"):
        check2(mut2(lambda s: s["modes"]["regime_cutoffs_ghz"].pop()))
    with pytest.raises(ValueError, match=r"^BBR025: .*regime_cutoffs_ghz"):
        check2(mut2(lambda s: s["modes"]["regime_cutoffs_ghz"].append({"mode": "TE07", "cutoff_ghz": 20983.0})))
    with pytest.raises(ValueError, match=r"^BBR025: .*TE0,2"):
        check2(mut2(lambda s: setp(s, "modes/regime_cutoffs_ghz/2/mode", "TE0,2")))
    # mistyped or absent: BBR024
    with pytest.raises(ValueError, match=r"^BBR024: missing field modes\.regime_cutoffs_ghz"):
        check2(mut2(lambda s: s["modes"].pop("regime_cutoffs_ghz")))
    with pytest.raises(ValueError, match=r"^BBR024: .*modes\.regime_cutoffs_ghz must be an array"):
        check2(mut2(lambda s: s["modes"].update(regime_cutoffs_ghz="TE10")))
    with pytest.raises(ValueError, match=r"^BBR024: .*mode"):
        check2(mut2(lambda s: s["modes"]["regime_cutoffs_ghz"][0].pop("mode")))
    # a section taller along g than l: TE01 is the lowest, so the ascending list starts with it
    tall = fx.sidecar2(y_half=1e-5, z_half=2e-5)
    check2(tall)
    cuts = [e["cutoff_ghz"] for e in tall["modes"]["regime_cutoffs_ghz"]]
    assert tall["modes"]["regime_cutoffs_ghz"][0]["mode"] == "TE01" and cuts == sorted(cuts)
    # the dropped schema-1 mode fields are not needed; carrying them is harmless (unknown fields ignored)
    def legacy(s):
        s["modes"].update(propagating_count=33, basis="closed PEC", list_limit_ghz=20000.0, mode_count_below_limit=9999,
                          gap_family_onsets=[{"n": 1, "mode": "TE01", "cutoff_ghz": 1.0}],
                          cutoffs=[{"mode": "TE11", "cutoff_ghz": 1.0, "degeneracy": 2}])
    check2(mut2(legacy))
    assert sidecar.same_invariant(mut2(legacy), fx.sidecar2())


def test_s2_cross_term_limit():
    lim = fx.sidecar2()["modes"]["polarization_filter_limit_ghz"]
    assert fx.sidecar2()["cross_term"]["valid_below_ghz"] == lim and lim == pytest.approx(2997.92458, rel=1e-6)
    with pytest.raises(ValueError, match=r"^BBR025: .*cross_term\.valid_below_ghz"):
        check2(mut2(lambda s: setp(s, "cross_term/valid_below_ghz", lim * (1 + 2e-6))))
    with pytest.raises(ValueError, match=r"^BBR025: .*polarization_filter_limit_ghz"):
        check2(mut2(lambda s: setp(s, "cross_term/valid_below_ghz", lim * (1 - 2e-6))))
    check2(mut2(lambda s: setp(s, "cross_term/valid_below_ghz", lim * (1 + 5e-7))))
    with pytest.raises(ValueError, match=r"^BBR024: .*cross_term\.valid_below_ghz"):
        check2(mut2(lambda s: s["cross_term"].pop("valid_below_ghz")))


def test_invariant_schema_major_and_grid():
    # Choice 1: one ID must not mix schema majors or incidence grids across frequencies
    one = built_1x()
    two = fx.sidecar2()
    d = sidecar.invariant_diff(one, two)
    assert d and d[0] == "schema" and not sidecar.same_invariant(one, two) and not sidecar.same_invariant(two, one)
    assert sidecar.invariant_diff(two, one) == sidecar.invariant_diff(one, two)
    assert sidecar.invariant(two)["schema"] == 2 and sidecar.invariant(one)["schema"] == 1
    # a 15-key style and a 35-key style grid of one geometry differ in "excitation" only
    a = fx.sidecar2()
    b = fx.sidecar2(thetas=(90, 105, 120, 135, 150, 165, 180), phis=(0, 22.5, 45, 67.5, 90))
    assert sidecar.invariant_diff(a, b) == ["excitation"] == sidecar.invariant_diff(b, a)
    assert not sidecar.same_invariant(a, b)
    # the same grid at another frequency agrees (2.x)
    assert sidecar.invariant_diff(a, fx.sidecar2()) == [] and sidecar.same_invariant(a, fx.sidecar2())
    # a changed phi list alone is also a different grid
    assert sidecar.invariant_diff(a, fx.sidecar2(phis=(0, 90))) == ["excitation"]
    # names come in order: schema first, then the existing blocks
    sec = fx.sidecar2(y_half=4.0e-3)
    assert sidecar.invariant_diff(one, sec)[0] == "schema"
    # the section also moves the extent_mm l that the check ties to it (CDX-028)
    assert sidecar.invariant_diff(a, sec) == ["geometry", "modes", "cross_section"]


def built_1x():
    """A schema-1 sidecar whose grid is the one of the schema-2 fixture's keys (0, 180) only."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        return tiny_sidecar(d, write_tiny(d))


# --- the 2.x CSV checks ---------------------------------------------------------------------------

def test_s2_csvs_pass(tmp_path):
    base, stem = s2_dataset(tmp_path)
    sc = fx.sidecar2()
    sidecar.check_csvs(base, stem, sc)
    base2 = str(tmp_path / "disc")
    os.makedirs(base2)
    stem2 = fx.write_nine_key_grid2(base2, disc=True)
    sidecar.check_csvs(base2, stem2, fx.sidecar2(disc=True))


def test_s2_header_bbr013(tmp_path):
    base, stem = s2_dataset(tmp_path)
    sc = fx.sidecar2()
    p = wg_path(base, stem, 0)
    text = open(p).read()
    # the 1.0 labels over 2.0 data: the header differs from the declared 2.0 list
    open(p, "w").write(text.replace(",".join(sidecar.WG_COLUMNS_2), ",".join(sidecar.WG_COLUMNS_2[:9] + [
        "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag", "InterceptFactor"]), 1))
    with pytest.raises(ValueError, match=r"^BBR013"):
        sidecar.check_csvs(base, stem, sc)
    # a 15-column table (no InterceptFactor) under a 2.x sidecar
    base2 = str(tmp_path / "b")
    os.makedirs(base2)
    fx.write_nine_key_grid2(base2)
    p = wg_path(base2, stem, 1)
    lines = open(p).read().splitlines()
    open(p, "w").write("\n".join(",".join(l.split(",")[:15]) for l in lines) + "\n")
    with pytest.raises(ValueError, match=r"^BBR013"):
        sidecar.check_csvs(base2, stem, sc)
    # a row with a field too many
    base3 = str(tmp_path / "c")
    os.makedirs(base3)
    fx.write_nine_key_grid2(base3)
    p = wg_path(base3, stem, 0)
    lines = open(p).read().splitlines()
    lines[3] += ",1"
    open(p, "w").write("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match=r"^BBR013"):
        sidecar.check_csvs(base3, stem, sc)


def test_s2_intercept_column_bbr012(tmp_path):
    sc = fx.sidecar2()
    scale = lambda f: (lambda v, n: v.__setitem__(15, repr(float(v[15]) * f)))
    # one row (the fifth) of (45, 135) Ephi=0 off by 2e-6: not constant within the key
    base, stem = s2_dataset(tmp_path / "a")
    edit_key(wg_path(base, stem, 0), 45, 135, lambda v, n: scale(1 + 2e-6)(v, n) if n == 4 else None)
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(45, 135\) polarization 0.*not constant within the key"):
        sidecar.check_csvs(base, stem, sc)
    # every row of the key off by 2e-6: constant, but not the sidecar's c_key
    base, stem = s2_dataset(tmp_path / "b")
    edit_key(wg_path(base, stem, 1), 45, 135, scale(1 + 2e-6))
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(45, 135\) polarization 1.*transmittance\.c_key"):
        sidecar.check_csvs(base, stem, sc)
    # within 1e-6 (5e-7, on one row and on every row of both polarizations): passes
    base, stem = s2_dataset(tmp_path / "c")
    edit_key(wg_path(base, stem, 0), 45, 135, lambda v, n: scale(1 + 5e-7)(v, n) if n == 4 else None)
    for e in (0, 1):
        edit_key(wg_path(base, stem, e), 90, 90, scale(1 + 5e-7))
    sidecar.check_csvs(base, stem, sc)


def test_s2_intercept_column_closed_form(tmp_path):
    """A sidecar whose c_key was edited alongside the column still has the reader's own recomputation to meet:
    the sidecar itself is refused (BBR025); a CSV column at the closed form's 2e-6 off while c_key is
    the one the sidecar declares is the c_key branch above. Both stand."""
    s = mut2(lambda s: setp(s, "transmittance/c_key/1/c_k", math.cos(math.radians(45)) * (1 + 2e-6)))
    with pytest.raises(ValueError, match=r"^BBR025"):
        check2(s)


def test_s2_two_t_sources_bbr012(tmp_path):
    # t_flux 0.6 in the sidecar, OutgoingPower 0.61e-9 over IngoingPower 1e-9 in the CSV
    base, stem = s2_dataset(tmp_path)
    edit_key(wg_path(base, stem, 0), 90, 180, lambda v, n: v.__setitem__(4, "6.1e-10"))
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(90, 180\) polarization 0.*t_flux / c_k"):
        sidecar.check_csvs(base, stem, fx.sidecar2())


def shift_key(base, stem, phi, theta, dphi):
    """Move the IWavePhi of key (phi, theta) by dphi in every CSV of both polarizations (a producer's key
    that is off its declared angle: within C2's 0.01 degree rounding, outside the sidecar's 1e-9)."""
    for e in (0, 1):
        for name in ("far_field.csv", "waveguide.csv"):
            edit_key(os.path.join(base, f"{stem}_Ephi={e}", name), phi, theta,
                     lambda v, n: v.__setitem__(2, repr(float(v[2]) + dphi)))


def test_s2_key_without_record_bbr012(tmp_path):
    """C2 matches keys after rounding to 0.01 degree; the per-key record is matched at 1e-9 degree. A CSV key
    in between has no record: refused (BBR012), not a silent skip of the row checks and the dark masking."""
    base, stem = s2_dataset(tmp_path / "a")
    shift_key(base, stem, 45, 135, 1e-7)
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(45\.0000001, 135\) polarization 0.*no per-key record"):
        sidecar.check_csvs(base, stem, fx.sidecar2())
    # the dark key's angle off by 1e-7 would otherwise lose the OR-038 check on its rows
    base, stem = s2_dataset(tmp_path / "b")
    shift_key(base, stem, 0, 180, 1e-7)
    with pytest.raises(ValueError, match=r"^BBR012: .*no per-key record"):
        sidecar.check_csvs(base, stem, fx.sidecar2())
    # 1e-10 degree is inside the sidecar's rule: passes
    base, stem = s2_dataset(tmp_path / "c")
    shift_key(base, stem, 45, 135, 1e-10)
    sidecar.check_csvs(base, stem, fx.sidecar2())


def test_s2_negative_outgoing_bbr012(tmp_path):
    base, stem = s2_dataset(tmp_path)
    edit_key(wg_path(base, stem, 0), 45, 135, lambda v, n: v.__setitem__(4, "-1e-12") if n == 3 else None)
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(45, 135\) polarization 0.*OutgoingPower -1e-12 is negative"):
        sidecar.check_csvs(base, stem, fx.sidecar2())


def test_s2_negative_tflux_bbr024():
    with pytest.raises(ValueError, match=r"^BBR024: .*polarization 0.*t_flux -0\.1 is negative"):
        check2(mut2(lambda s: setp(s, "transmittance/per_key/1/polarizations/0/t_flux", -0.1)))
    with pytest.raises(ValueError, match=r"^BBR024: .*polarization 1.*is negative"):
        check2(mut2(lambda s: setp(s, "transmittance/per_key/3/polarizations/1/t_flux", -1e-300)))


def test_s2_dark_outgoing_bbr012(tmp_path):
    # OR-038: a dark basis needs OutgoingPower exactly 0 on every row; 1e-25 on the third row is refused
    base, stem = s2_dataset(tmp_path)
    edit_key(wg_path(base, stem, 1), 0, 180, lambda v, n: v.__setitem__(4, "1e-25") if n == 2 else None)
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(0, 180\) polarization 1.*OutgoingPower = (1e-25|9\.99\d*e-26)"):
        sidecar.check_csvs(base, stem, fx.sidecar2())
    # the dark basis keeps nonzero stored fields (the fixture's), which is fine
    base2 = str(tmp_path / "ok")
    os.makedirs(base2)
    fx.write_nine_key_grid2(base2)
    sidecar.check_csvs(base2, stem, fx.sidecar2())


def test_s2_dark_first_row_residue_bbr012(tmp_path):
    """A dark basis whose first row is the residue (T = 1e-16 here) contradicts the sidecar's t_flux = 0:
    both T sources must be exactly 0 on a dark key."""
    base, stem = s2_dataset(tmp_path)
    edit_key(wg_path(base, stem, 1), 0, 180, lambda v, n: v.__setitem__(4, "1e-25"))
    with pytest.raises(ValueError, match=r"^BBR012: .*key \(0, 180\) polarization 1"):
        sidecar.check_csvs(base, stem, fx.sidecar2())


# --- the schema-2 builder -------------------------------------------------------------------------

def build2(base, stem, **kw):
    ds_id, label = stem.rsplit("_", 1)
    section = kw.pop("section", {"shape": "rectangle", "y_e_half_m": 4.5e-3, "z_e_half_m": 2.5e-5})
    width = ((2e3 * section["radius_m"],) * 2 if section["shape"] == "disc"
             else (2e3 * section["y_e_half_m"], 2e3 * section["z_e_half_m"]))
    return sidecar.build_schema2_from_csvs(
        base, ds_id, label, section=section, extent_mm={"p": 4.0, "l": width[0], "g": width[1]},
        provenance={"producer": "test", "material": "vacuum"}, exit_origin_mm=[0.0, 0.0, 4.0],
        plane_wave_origin_mm=[0.0, 0.0, 0.0], bounding_box_mm=[-0.026, -5.0, 0.0, 0.026, 5.0, 4.0], **kw)


def digest(base):
    out = {}
    for root, _, files in os.walk(base):
        for f in files:
            with open(os.path.join(root, f), "rb") as fh:
                out[os.path.join(root, f)] = hashlib.sha256(fh.read()).hexdigest()
    return out


def test_s2_builder_roundtrip(tmp_path):
    base, stem = s2_dataset(tmp_path)
    before = digest(base)
    sc = build2(base, stem)
    assert digest(base) == before                                     # the builder never modifies its inputs
    sidecar.write(base, stem, sc)
    got = sidecar.check_full(base, ID, stem, 500.0)
    assert got["schema_version"] == "2.0" and sidecar.schema_major(got) == 2
    ref = fx.sidecar2()
    # the blocks the CSVs and the declared grid determine equal the independent reference
    for blk in ("symmetry", "boundaries", "model", "conventions", "dataset_id", "frequency_label", "frequency_ghz"):
        assert got[blk] == ref[blk], blk
    # the two copies of the closed forms agree to rounding (the fixture's are written independently)
    assert got["modes"]["mode"] == ref["modes"]["mode"]
    assert got["modes"]["cutoff_ghz"] == pytest.approx(ref["modes"]["cutoff_ghz"], rel=1e-12)
    assert got["modes"]["polarization_filter_limit_ghz"] == pytest.approx(
        ref["modes"]["polarization_filter_limit_ghz"], rel=1e-12)
    assert [e["mode"] for e in got["modes"]["regime_cutoffs_ghz"]] == [e["mode"] for e in ref["modes"]["regime_cutoffs_ghz"]]
    assert [e["cutoff_ghz"] for e in got["modes"]["regime_cutoffs_ghz"]] == pytest.approx(
        [e["cutoff_ghz"] for e in ref["modes"]["regime_cutoffs_ghz"]], rel=1e-12)
    assert {k: v for k, v in got["cross_term"].items() if k != "valid_below_ghz"} == \
        {k: v for k, v in ref["cross_term"].items() if k != "valid_below_ghz"}
    assert got["cross_term"]["valid_below_ghz"] == pytest.approx(ref["cross_term"]["valid_below_ghz"], rel=1e-12)
    assert {k: got["excitation"][k] for k in ("incident_phi_deg", "incident_theta_deg")} == \
        {k: ref["excitation"][k] for k in ("incident_phi_deg", "incident_theta_deg")}
    assert got["frames"]["hfss_global_axes_in_canonical"] == ref["frames"]["hfss_global_axes_in_canonical"]
    t, rt = got["transmittance"], ref["transmittance"]
    assert t["c_key"] == rt["c_key"] and t["cell_90"] == rt["cell_90"] and t["cell_90_assumes"] == rt["cell_90_assumes"]
    for g, r in zip(t["per_key"], rt["per_key"]):
        assert (g["phi_deg"], g["theta_deg"], g["c_k"]) == (r["phi_deg"], r["theta_deg"], r["c_k"])
        for p in "01":
            assert g["polarizations"][p]["t_flux"] == pytest.approx(r["polarizations"][p]["t_flux"], rel=1e-12)
    # no other numerator given: the magnitude is the stored OutgoingPower
    assert t["per_key"][0]["polarizations"]["0"]["outgoing_power_abs_re_s_w"] == pytest.approx(0.6e-9)
    scale = got["exit_field"]["scale"]["per_key"]
    assert all(s["w0"] == 1.0 and s["w1"] == 1.0 for s in scale)
    dark = [s for s in scale if s["dark1"]]
    assert len(dark) == 1 and (dark[0]["phi_deg"], dark[0]["theta_deg"]) == (0, 180) and not dark[0]["dark0"]
    ffr = got["far_field"]["ff_over_flux"]
    for rec in ffr:
        for p in "01":
            v = rec["polarizations"][p]
            dk = (rec["phi_deg"], rec["theta_deg"], p) == (0, 180, "1")
            assert (v is None) == dk and (dk or (math.isfinite(v) and v > 0))
    assert got["exit_field"]["grid"]["y_e"]["step"] == pytest.approx(1e-3)
    assert got["exit_field"]["grid"]["z_e"]["step"] == pytest.approx(2.5e-5)
    assert got["exit_field"]["columns"] == sidecar.WG_COLUMNS_2 and got["exit_field"]["quadrature"] == "trapezoid"
    assert got["files"] and all(k.startswith(stem) for k in got["files"])


def test_s2_builder_w_and_steps(tmp_path):
    base, stem = s2_dataset(tmp_path)
    w = {(ph, t): (0.9, 0.8) for ph in fx.PHIS for t in fx.THETAS}
    sc = build2(base, stem, w=w, step_m={"y_e": 1e-3, "z_e": 2.5e-5}, rotational=True)
    scale = {(s["phi_deg"], s["theta_deg"]): s for s in sc["exit_field"]["scale"]["per_key"]}
    assert (scale[(45, 135)]["w0"], scale[(45, 135)]["w1"]) == (0.9, 0.8)
    assert (scale[(0, 180)]["w0"], scale[(0, 180)]["w1"]) == (0.9, 1.0)       # a dark basis has W = 1
    assert sc["symmetry"]["rotational"] is True
    sidecar.write(base, stem, sc)
    sidecar.check_full(base, ID, stem, 500.0)
    # a scalar step applies to both axes (the step must still meet the lattice, so a wrong one is refused)
    sc2 = build2(base, stem, step_m=1e-3)
    with pytest.raises(ValueError, match=r"^BBR024: .*step"):
        sidecar.check(sc2, ID, stem, 500.0)
    # abs_powers override the copied numerator
    sc3 = build2(base, stem, abs_powers={(0, 180): (7e-10, 5e-22)})
    pk = [p for p in sc3["transmittance"]["per_key"] if (p["phi_deg"], p["theta_deg"]) == (0, 180)][0]
    assert pk["polarizations"]["0"]["outgoing_power_abs_re_s_w"] == 7e-10
    assert pk["polarizations"]["1"]["outgoing_power_abs_re_s_w"] == 5e-22


def test_s2_builder_disc_roundtrip(tmp_path):
    base, stem = s2_dataset(tmp_path, disc=True)
    sc = build2(base, stem, section={"shape": "disc", "radius_m": 5e-5})
    sidecar.write(base, stem, sc)
    got = sidecar.check_full(base, ID, stem, 500.0)
    assert got["exit_field"]["quadrature"] == "uniform_disc" and got["modes"]["mode"] == "TE11"
    assert [e["mode"] for e in got["modes"]["regime_cutoffs_ghz"]] == ["TE11", "TM01", "TE21"]
    assert got["geometry"]["shape"] == "cylinder"
    assert got["cross_term"]["valid_below_ghz"] == got["modes"]["polarization_filter_limit_ghz"]
    # and it carries no 1.x mode lists
    assert not {"gap_family_onsets", "cutoffs", "propagating_count", "basis", "list_limit_ghz"} & set(got["modes"])


def test_s2_builder_refuses_a_dark_basis_with_residue(tmp_path):
    base, stem = s2_dataset(tmp_path)
    # one row of the dark basis holds a numerical residue: refused, nothing written
    edit_key(wg_path(base, stem, 1), 0, 180, lambda v, n: v.__setitem__(4, "1e-25") if n == 2 else None)
    with pytest.raises(ValueError, match=r"^BBR012: .*dark.*OutgoingPower"):
        build2(base, stem)
    # the first row alone being the residue is also refused (|t_flux| = 1e-16 <= 1e-12 but nonzero)
    base2 = str(tmp_path / "b")
    os.makedirs(base2)
    fx.write_nine_key_grid2(base2)
    edit_key(wg_path(base2, stem, 1), 0, 180, lambda v, n: v.__setitem__(4, "1e-25"))
    with pytest.raises(ValueError, match=r"^BBR012: .*dark.*OutgoingPower"):
        build2(base2, stem)
    assert not os.path.exists(sidecar.path_for(base2, stem))


def test_s2_builder_needs_the_2x_header(tmp_path):
    base = str(tmp_path)
    stem = write_tiny(base)                                           # 1.0 tables
    with pytest.raises(ValueError, match=r"^BBR013"):
        build2(base, stem, section={"shape": "rectangle", "y_e_half_m": 2e-3, "z_e_half_m": 2.5e-5})


def test_s2_check_full_on_a_stale_checksum(tmp_path):
    base, stem = s2_dataset(tmp_path)
    sc = build2(base, stem)
    sidecar.write(base, stem, sc)
    with open(wg_path(base, stem, 0), "a") as fh:
        fh.write("\n")
    with pytest.raises(ValueError, match=r"^BBR024: .*sha256"):
        sidecar.check_full(base, ID, stem, 500.0)
