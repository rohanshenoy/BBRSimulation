"""End-to-end checks of the per-track reflection validator on ROOT output."""
import json
import os
import subprocess
import sys
import importlib.util

import numpy as np
import pandas as pd
import pytest
import uproot

from bbrsim.io import load


LEGEND = {
    "status": {"0": "FresnelRefraction", "1": "FresnelReflection",
               "2": "TIR", "14": "BBRDiffractionTransmit",
               "15": "BBRDiffractionReflect", "16": "BBRReflect",
               "12": "NoRINDEX", "13": "Other", "17": "BBRAbsorb", "18": "unknown", "19": "WorldExit",
               "21": "CoatedDielectricReflection", "22": "CoatedDielectricRefraction",
               "23": "CoatedDielectricFrustratedTransmission",
               "24": "PolishedLumirrorAirReflection"},
    "event_type": {"0": "transmission", "1": "reflection", "2": "absorption", "3": "other"},
    "volume": {"-1": "none", "0": "World"},
    "material": {"-1": "none", "0": "G4_Galactic"},
}


def crossing(event_id, boundary, reflections, status, track_id=1, run_id=0):
    return dict(run_id=run_id, event_id=event_id, track_id=track_id,
                n_boundary=boundary, n_reflections=reflections,
                status_code=status)


def termination(event_id, boundary, reflections, status=19, track_id=1, run_id=0):
    return dict(run_id=run_id, event_id=event_id, track_id=track_id,
                n_boundary=boundary, n_reflections=reflections,
                term_status_code=status)


def write_output(directory, crossings, terminations, *, legacy=False):
    directory.mkdir(exist_ok=True)
    cr_names = ("run_id", "event_id", "track_id", "n_boundary", "n_reflections",
                "status_code", "event_type_code", "vol_pre_code", "vol_post_code",
                "mat_pre_code", "mat_post_code")
    ab_names = ("run_id", "event_id", "track_id", "n_boundary", "n_reflections",
                "term_status_code", "term_vol_code")
    if legacy:
        cr_names = tuple(n for n in cr_names if n not in ("track_id", "n_boundary", "n_reflections"))
        ab_names = tuple(n for n in ab_names if n not in ("track_id", "n_boundary", "n_reflections"))
    cr = {n: np.asarray([row.get(n, 0) for row in crossings], dtype=np.int32) for n in cr_names}
    ab = {n: np.asarray([row.get(n, 0) for row in terminations], dtype=np.int32) for n in ab_names}
    path = directory / "bbr.root"
    with uproot.recreate(path) as root:
        root.mktree("crossings", {n: arr.dtype for n, arr in cr.items()})
        if crossings:
            root["crossings"].extend(cr)
        root.mktree("abspoints", {n: arr.dtype for n, arr in ab.items()})
        if terminations:
            root["abspoints"].extend(ab)
    if legacy:
        (directory / "bbr_legend.json").write_text(json.dumps(LEGEND))
    else:
        path.with_suffix(".metadata.json").write_text(json.dumps({
            "schema_version": 2, "legend": LEGEND, "run_id": 0,
            "configuration": "synthetic reflection-count fixture",
            "build": {"version": "test"},
        }))
    return path


def check(repo_root, tmp_path, crossings, terminations, *, legacy=False):
    path = write_output(tmp_path / "out", crossings, terminations, legacy=legacy)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "tools/python") + os.pathsep + env.get("PYTHONPATH", "")
    env["MPLCONFIGDIR"] = str(tmp_path / "mpl")
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    return subprocess.run([sys.executable, str(repo_root / "validation/check_nreflect.py"), str(path)],
                          cwd=tmp_path, env=env, text=True, capture_output=True)


@pytest.fixture(scope="module")
def validator(repo_root):
    spec = importlib.util.spec_from_file_location("check_nreflect", repo_root / "validation/check_nreflect.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_out_of_order_crossings_and_transmission_are_valid(repo_root, tmp_path):
    # A transmission raises the crossing ordinal but does not add a reflection.
    cr = [crossing(0, 3, 1, 17), crossing(0, 1, 0, 14), crossing(0, 2, 1, 16),
          crossing(1, 1, 1, 15, track_id=2)]
    ab = [termination(1, 1, 1, track_id=2), termination(0, 3, 1, status=17)]
    result = check(repo_root, tmp_path, cr, ab)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "RESULT: PASS" in result.stdout
    assert (tmp_path / "nreflect_distribution.png").exists()


def test_zero_crossing_world_exit_is_valid(repo_root, tmp_path):
    result = check(repo_root, tmp_path, [], [termination(0, 0, 0)])
    assert result.returncode == 0, result.stdout + result.stderr
    assert "RESULT: PASS" in result.stdout


def test_same_event_and_track_ids_in_different_runs_are_distinct(repo_root, tmp_path):
    cr = [crossing(0, 1, 1, 16, run_id=0), crossing(0, 1, 0, 14, run_id=1)]
    ab = [termination(0, 1, 0, run_id=1), termination(0, 1, 1, run_id=0)]
    result = check(repo_root, tmp_path, cr, ab)
    assert result.returncode == 0, result.stdout + result.stderr


def test_same_event_with_different_tracks_is_distinct(repo_root, tmp_path):
    cr = [crossing(0, 1, 1, 16, track_id=1), crossing(0, 1, 0, 14, track_id=2)]
    ab = [termination(0, 1, 0, track_id=2), termination(0, 1, 1, track_id=1)]
    result = check(repo_root, tmp_path, cr, ab)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("status_code", [21, 24], ids=["coated", "lut"])
def test_stock_coated_and_lut_reflections_counted(validator, tmp_path, status_code):
    path = write_output(tmp_path / "out", [crossing(0, 1, 1, status_code)],
                        [termination(0, 1, 1)])
    counts, errors = validator.validate(*load(path))
    assert errors == []
    assert counts == [1]


@pytest.mark.parametrize("status_code", [22, 23], ids=["coated-refraction", "frustrated-transmission"])
def test_coated_transmissions_do_not_count_as_reflections(validator, tmp_path, status_code):
    path = write_output(tmp_path / "out", [crossing(0, 1, 0, status_code)],
                        [termination(0, 1, 0)])
    counts, errors = validator.validate(*load(path))
    assert errors == []
    assert counts == [0]


@pytest.mark.parametrize("status", [
    "PolishedLumirrorAirReflection", "PolishedLumirrorGlueReflection",
    "PolishedAirReflection", "PolishedTeflonAirReflection",
    "PolishedTiOAirReflection", "PolishedTyvekAirReflection",
    "PolishedVM2000AirReflection", "PolishedVM2000GlueReflection",
    "EtchedLumirrorAirReflection", "EtchedLumirrorGlueReflection",
    "EtchedAirReflection", "EtchedTeflonAirReflection",
    "EtchedTiOAirReflection", "EtchedTyvekAirReflection",
    "EtchedVM2000AirReflection", "EtchedVM2000GlueReflection",
    "GroundLumirrorAirReflection", "GroundLumirrorGlueReflection",
    "GroundAirReflection", "GroundTeflonAirReflection",
    "GroundTiOAirReflection", "GroundTyvekAirReflection",
    "GroundVM2000AirReflection", "GroundVM2000GlueReflection",
])
def test_every_lut_reflection_is_counted(validator, status):
    cr = pd.DataFrame([dict(run_id=0, event_id=0, track_id=1,
                            n_boundary=1, n_reflections=1, status=status)])
    ab = pd.DataFrame([dict(run_id=0, event_id=0, track_id=1,
                            n_boundary=1, n_reflections=1, term_status="WorldExit")])
    assert validator.validate(cr, ab) == ([1], [])


@pytest.mark.parametrize("crossings,terminations,signature", [
    ([crossing(0, 1, 0, 16)], [termination(0, 1, 0)], "reflection count"),
    ([crossing(0, 1, 1, 14)], [termination(0, 1, 1)], "reflection count"),
    ([crossing(0, 1, 1, 16)], [termination(0, 1, 0)], "termination"),
    ([crossing(0, 1, 0, 14), crossing(0, 1, 1, 16)], [termination(0, 1, 1)], "duplicate"),
    ([crossing(0, 1, 0, 14), crossing(0, 3, 1, 16)], [termination(0, 3, 1)], "gap"),
    ([crossing(0, 1, 0, 18)], [termination(0, 1, 0)], "unknown status"),
    ([crossing(0, 1, 0, 13)], [termination(0, 1, 0)], "unknown status"),
    ([crossing(0, 1, 0, 99)], [termination(0, 1, 0)], "unknown status"),
    ([crossing(0, 1, 0, 14)], [], "termination"),
    ([], [termination(0, 0, 0), termination(0, 0, 0)], "duplicate termination"),
    ([], [termination(0, 0, 0, status=17)], "missing killing boundary"),
    ([], [termination(0, 0, 0, status=12)], "missing killing boundary"),
    ([crossing(0, 1, 0, 14)], [termination(0, 1, 0, status=17)], "killing boundary"),
], ids=["missed-reflection", "transmission-counted", "termination-mismatch",
        "duplicate-ordinal", "missing-ordinal", "unknown-label", "other-label", "missing-legend-code",
        "missing-termination", "duplicate-termination", "missing-killing-boundary",
        "missing-no-rindex-boundary", "wrong-killing-status"])
def test_corrupt_counts_fail(validator, tmp_path, crossings, terminations, signature):
    path = write_output(tmp_path / "out", crossings, terminations)
    cr, ab = load(path)
    _, errors = validator.validate(cr, ab)
    assert errors
    assert signature in "\n".join(errors).lower()


def test_empty_output_fails(repo_root, tmp_path):
    result = check(repo_root, tmp_path, [], [])
    assert result.returncode != 0
    assert "empty output" in result.stdout.lower()


def test_legacy_schema_fails_explicitly(repo_root, tmp_path):
    result = check(repo_root, tmp_path, [], [termination(0, 0, 0)], legacy=True)
    assert result.returncode != 0
    assert "unsupported schema" in result.stdout.lower()
