"""bbrsim.io on ROOT files uproot writes in the test (mktree + extend), plus select and paths."""
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest
import uproot

from bbrsim import io, paths, select

LEGEND = {
    "status": {"14": "BBRDiffractionTransmit", "16": "BBRReflect", "17": "BBRAbsorb", "19": "WorldExit"},
    "event_type": {"0": "transmission", "1": "reflection", "2": "absorption", "3": "other"},
    "volume": {"0": "CuSlab", "1": "InfParallelPlate_crack1Rohan", "3": "World", "-1": "none"},
    "material": {"0": "Cu_RRR100_T4K", "1": "G4_Galactic", "2": "vacuum_wg", "-1": "none"},
}
I4, F8 = np.int32, np.float64
CR_INT = ["run_id", "event_id", "vol_pre_code", "mat_pre_code", "vol_post_code", "mat_post_code",
          "status_code", "event_type_code", "n_reflect"]
CR_DBL = ["x_mm", "y_mm", "z_mm", "energy_eV", "px_pre", "py_pre", "pz_pre", "px_post", "py_post",
          "pz_post", "theta_in_deg", "phi_in_deg", "hfss_freq_GHz"]
AB_INT = ["run_id", "event_id", "n_reflect", "term_vol_code", "term_status_code"]
AB_DBL = ["x_mm", "y_mm", "z_mm", "energy_eV", "px", "py", "pz"]


def write_output(d, crossings, abspoints, run_id=0):
    """Write a G4Analysis-shaped bbr.root + bbr_legend.json into directory d."""
    os.makedirs(d, exist_ok=True)
    n = len(crossings["status_code"])
    cr = {c: np.asarray(crossings.get(c, [0] * n), I4) for c in CR_INT}
    cr.update({c: np.asarray(crossings.get(c, [0.0] * n), F8) for c in CR_DBL})
    cr["run_id"] = np.full(n, run_id, I4)
    m = len(abspoints.get("term_status_code", []))
    ab = {c: np.asarray(abspoints.get(c, [0] * m), I4) for c in AB_INT}
    ab.update({c: np.asarray(abspoints.get(c, [0.0] * m), F8) for c in AB_DBL})
    path = os.path.join(d, "bbr.root")
    with uproot.recreate(path) as f:
        f.mktree("crossings", {k: v.dtype for k, v in cr.items()})
        if n:
            f["crossings"].extend(cr)
        f.mktree("abspoints", {k: v.dtype for k, v in ab.items()})
        if m:
            f["abspoints"].extend(ab)
    with open(os.path.join(d, "bbr_legend.json"), "w") as fh:
        json.dump(LEGEND, fh)
    return path


CROSS = {"event_id": [0, 0, 1, 2], "n_reflect": [1, 2, 1, 1],
         "vol_pre_code": [3, 3, 3, 3], "mat_pre_code": [1, 1, 1, 1],
         "vol_post_code": [0, 0, 1, 0], "mat_post_code": [0, 0, 2, 0],
         "status_code": [16, 17, 14, -2], "event_type_code": [1, 2, 0, 3]}
ABS = {"event_id": [0, 1], "term_vol_code": [0, -1], "term_status_code": [17, 19]}


def test_load_decodes_codes(tmp_path):
    p = write_output(str(tmp_path / "out"), CROSS, ABS)
    cr, ab = io.load(p)
    assert cr["status"].tolist()[:3] == ["BBRReflect", "BBRAbsorb", "BBRDiffractionTransmit"]
    assert pd.isna(cr["status"].iloc[3])            # code -2 (missing legend entry) -> NaN
    assert cr["vol_post"].tolist() == ["CuSlab", "CuSlab", "InfParallelPlate_crack1Rohan", "CuSlab"]
    assert ab["term_vol"].tolist() == ["CuSlab", "none"] and ab["term_status"].tolist() == ["BBRAbsorb", "WorldExit"]
    assert cr["run_id"].dtype == np.int32 and cr["x_mm"].dtype == np.float64


def test_load_empty_abspoints_has_empty_decoded_columns(tmp_path):
    # D18: the decoded term_* columns exist even when abspoints has no rows.
    p = write_output(str(tmp_path / "out"), CROSS, {})
    _, ab = io.load(p)
    assert len(ab) == 0 and {"term_status", "term_vol"} <= set(ab.columns)


def test_load_missing_legend_raises(tmp_path):
    p = write_output(str(tmp_path / "out"), CROSS, ABS)
    os.remove(os.path.join(os.path.dirname(p), "bbr_legend.json"))
    with pytest.raises(FileNotFoundError):
        io.load(p)


def test_load_many_keeps_runs_apart(tmp_path):
    p0 = write_output(str(tmp_path / "r0"), CROSS, ABS, run_id=0)
    p1 = write_output(str(tmp_path / "r1"), CROSS, ABS, run_id=1)
    df = io.load_many([p0, p1])
    assert len(df) == 8 and sorted(df["run_id"].unique()) == [0, 1]
    assert io.load_many([]).empty


def test_select_helpers(tmp_path):
    df = io.load_crossings(write_output(str(tmp_path / "out"), CROSS, ABS))
    assert len(select.cu_boundary(df)) == 3
    assert len(select.first_hit_cu(df)) == 2
    assert select.crack_crossings(df)["vol_post"].tolist() == ["InfParallelPlate_crack1Rohan"]
    assert select.cu_absorption_stats(df) == (2, 1)
    assert "evt_key" not in df.columns                 # cu_absorption_stats copies
    assert select.add_evt_key(df)["evt_key"].iloc[0] == (0, 0)


@pytest.mark.parametrize("name,expected", [
    ("Cu_RRR100_T4K", (100, 4.0)), ("Cu_RRR50_T77.5K", (50, 77.5)), ("G4_Galactic", None)])
def test_parse_cu_rrr_t(name, expected):
    assert select.parse_cu_rrr_t(name) == expected


def test_paths_env_wins_and_empty_counts_as_unset(monkeypatch, tmp_path):
    monkeypatch.setenv("BBRSIMDATA", "/nonexistent/x")
    assert paths.data_dir() == "/nonexistent/x"
    fake = tmp_path / "tree" / "tools" / "python" / "bbrsim"
    fake.mkdir(parents=True)
    (tmp_path / "tree" / "data" / "waveguides").mkdir(parents=True)
    monkeypatch.setattr(paths, "__file__", str(fake / "paths.py"))
    monkeypatch.setenv("BBRSIMDATA", "")
    assert paths.data_dir() == str(tmp_path / "tree" / "data")


def test_paths_prefix_fallback_and_error(monkeypatch, tmp_path):
    monkeypatch.delenv("BBRSIMDATA", raising=False)
    monkeypatch.setattr(paths, "__file__", str(tmp_path / "a" / "paths.py"))
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "pfx"))
    with pytest.raises(FileNotFoundError):
        paths.data_dir()
    (tmp_path / "pfx" / "share" / "BBRsim" / "data").mkdir(parents=True)
    assert paths.data_dir() == str(tmp_path / "pfx" / "share" / "BBRsim" / "data")
