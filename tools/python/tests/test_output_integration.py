"""Real executable checks; enabled by the regression runner's binary paths."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from bbrsim.io import load
from bbrsim.paths import data_dir


@pytest.fixture
def binary():
    name = os.environ.get("BBR_TESTWORLD_BINARY")
    if not name:
        pytest.skip("test-world binary not supplied")
    return name


@pytest.fixture(scope="module")
def round_gap_data(repo_root, tmp_path_factory):
    """A data root holding the mock round gap beside links to the real cracks."""
    if not os.environ.get("BBR_TESTWORLD_BINARY"):   # its tests skip; do not build it
        pytest.skip("test-world binary not supplied")
    root = tmp_path_factory.mktemp("round_gap_data")
    gen = repo_root / "validation/Scripts/make_mock_round_gap.py"
    r = subprocess.run([sys.executable, str(gen), "--real", str(Path(data_dir()) / "waveguides"),
                        "--dst", str(root)], capture_output=True, text=True, cwd=root)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (root / "waveguides" / "RoundGap_r50um_2000GHz.dataset.json").is_file()
    return root


def run_macro(binary, work, text):
    macro = work / "test.mac"
    macro.write_text(text)
    env = dict(os.environ, G4FORCENUMBEROFTHREADS="1")
    env.pop("DYLD_LIBRARY_PATH", None)
    env.pop("LD_LIBRARY_PATH", None)
    return subprocess.run([binary, str(macro)], cwd=work, env=env,
                          capture_output=True, text=True, timeout=60)


def commands(name, extra=""):
    return (extra + "\n/run/initialize\n/bbr/gun/mode true\n/bbr/gun/posY 1\n"
            "/bbr/gun/posZ 1\n/bbr/gun/dirX 1\n"
            f"/analysis/setFileName {name}\n/run/beamOn 20\n")


def test_redirected_results_keep_independent_metadata(binary, tmp_path):
    first = run_macro(binary, tmp_path, commands("results/testworld.root"))
    assert first.returncode == 0, first.stdout + first.stderr
    root = tmp_path / "results/testworld.root"
    root_bytes = root.read_bytes()
    meta_bytes = root.with_suffix(".metadata.json").read_bytes()
    cr_before, ab_before = load(root)
    # A second, crack-free run of the same binary in the same directory.
    second = run_macro(binary, tmp_path, commands("results/nocracks.root",
                                                  "/bbr/testworld/cracks false"))
    assert second.returncode == 0, second.stdout + second.stderr
    assert root.read_bytes() == root_bytes
    assert root.with_suffix(".metadata.json").read_bytes() == meta_bytes
    cr_after, ab_after = load(root)
    assert cr_before.equals(cr_after) and ab_before.equals(ab_after)
    assert len(ab_after) == 20
    meta = json.loads(root.with_suffix(".metadata.json").read_text())
    assert meta["schema_version"] == 2 and meta["run_id"] == 0
    assert "/bbr/gun/mode" in meta["configuration"]
    assert meta["build"]["version"] and meta["build"]["source_fingerprint"]
    assert meta["build"]["application_source_fingerprint"].startswith("sha256:")
    assert meta["data"]["directory"] and "fingerprint" in meta["data"]
    assert {"track_id", "n_boundary", "n_reflections"} <= set(cr_after)
    assert {"track_id", "n_boundary", "n_reflections"} <= set(ab_after)
    assert not cr_after["vol_post"].isna().any()
    # The two placed cracks, their extents and each frequency's recorded sidecar fields.
    cracks = {c["dataset_id"]: c for c in meta["data"]["hfss_datasets"]}
    assert len(meta["data"]["hfss_datasets"]) == 2
    assert set(cracks) == {"InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2"}
    expected = {"InfParallelPlate_crack1Rohan": ([1.0, 10.0, 0.05], {"p": 1.0, "l": 10.0, "g": 0.05}),
                "InfParallelPlate_crack2": ([1.5, 10.0, 0.1], {"p": 1.5, "l": 10.0, "g": 0.1})}
    for dataset_id, (geant4, hfss) in expected.items():
        crack = cracks[dataset_id]
        assert crack["geant4_extent_mm"] == pytest.approx(geant4, abs=1e-9)
        assert crack["hfss_extent_mm"] == hfss
        assert len(crack["frequencies"]) == 1
        frequency = crack["frequencies"][0]
        assert frequency["label"] == "500GHz"
        assert frequency["sidecar"] == f"{dataset_id}_500GHz.dataset.json"
        assert frequency["recorded"]["provenance"]["producer"] == "legacy bbsim1freq.py"
        assert frequency["recorded"]["frames"]["pose_rule"] == "legacy"
    nocracks = json.loads((tmp_path / "results/nocracks.metadata.json").read_text())
    assert nocracks["data"]["hfss_datasets"] == []
    assert "CuPlateCrack" not in nocracks["geometry"]   # the plates go with their cracks


def test_output_failure_is_fatal(binary, tmp_path):
    (tmp_path / "blocked").write_text("a file, not a directory")
    result = run_macro(binary, tmp_path, commands("blocked/test.root"))
    assert result.returncode != 0, result.stdout + result.stderr
    assert "BBR022" in result.stdout + result.stderr


def test_multirun_metadata_captures_configuration_and_json_escaping(binary, tmp_path):
    data = tmp_path / 'data"quote'
    data.mkdir()
    # The test world's cracks are validated at /run/initialize, so the root needs their datasets.
    (data / "waveguides").symlink_to(Path(data_dir()) / "waveguides", target_is_directory=True)
    macro = (f"/bbr/dataDir '{data}'\n" + commands("one.root") +
             "/bbr/gun/energy_eV 0.003\n/analysis/setFileName two.root\n/run/beamOn 10\n")
    result = run_macro(binary, tmp_path, macro)
    assert result.returncode == 0, result.stdout + result.stderr
    first = json.loads((tmp_path / "one.metadata.json").read_text())
    second = json.loads((tmp_path / "two.metadata.json").read_text())
    assert first["run_id"] == 0 and second["run_id"] == 1
    assert first["events"] == 20 and second["events"] == 10
    assert first["worker_count"] == 1
    assert first["configuration"] != second["configuration"]
    assert first["data"]["directory"] == str(data)
    assert first["geometry"] and first["random_engine_state"]
    # The crack records are read once per process; a later run must carry the same ones.
    assert len(first["data"]["hfss_datasets"]) == 2
    assert second["data"]["hfss_datasets"] == first["data"]["hfss_datasets"]


def test_round_gap_switch_is_refused_after_initialize(binary, tmp_path):
    # The round gap changes the geometry, which is built at /run/initialize: the
    # switch is PreInit only, and a batch run fails on the refused command.
    result = run_macro(binary, tmp_path, "/run/initialize\n/bbr/testworld/roundGap true\n")
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "Illegal application state </bbr/testworld/roundGap true>" in output


def test_cracks_switch_is_refused_after_initialize(binary, tmp_path):
    # The cracks are geometry, built at /run/initialize: the switch is PreInit only.
    result = run_macro(binary, tmp_path, "/run/initialize\n/bbr/testworld/cracks false\n")
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "Illegal application state </bbr/testworld/cracks false>" in output


def test_cracks_off_keeps_round_gap(binary, tmp_path, round_gap_data):
    # Only the two slab plates and their cracks go; the round gap is its own switch.
    macro = (f"/bbr/dataDir {round_gap_data}\n/bbr/testworld/cracks false\n"
             "/bbr/testworld/roundGap true\n" + commands("rg.root").replace("/run/beamOn 20", "/run/beamOn 5"))
    result = run_macro(binary, tmp_path, macro)
    assert result.returncode == 0, result.stdout + result.stderr
    meta = json.loads((tmp_path / "rg.metadata.json").read_text())
    assert [d["dataset_id"] for d in meta["data"]["hfss_datasets"]] == ["RoundGap_r50um"]
    assert "CuPlateCrack" not in meta["geometry"] and "RoundGapPlate" in meta["geometry"]


def test_cracks_off_needs_no_waveguides(binary, tmp_path):
    # No vacuum_wg volume, so nothing is discovered: a data root without waveguides/ runs.
    empty = tmp_path / "empty"
    empty.mkdir()
    macro = (f"/bbr/dataDir {empty}\n/bbr/testworld/cracks false\n"
             + commands("free.root").replace("/run/beamOn 20", "/run/beamOn 5"))
    result = run_macro(binary, tmp_path, macro)
    assert result.returncode == 0, result.stdout + result.stderr
    meta = json.loads((tmp_path / "free.metadata.json").read_text())
    assert meta["data"]["hfss_datasets"] == []
    assert "CuPlateCrack" not in meta["geometry"]
