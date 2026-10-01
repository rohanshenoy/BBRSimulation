"""Real executable checks; enabled by the regression runner's binary paths."""
import json
import os
from pathlib import Path
import subprocess

import pytest
from bbrsim.io import load


@pytest.fixture
def binaries():
    names = [os.environ.get("BBR_TESTWORLD_BINARY"), os.environ.get("BBR_LIGHTPIPE_BINARY")]
    if not all(names):
        pytest.skip("example binaries not supplied")
    return names


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


def test_redirected_results_keep_independent_metadata(binaries, tmp_path):
    first = run_macro(binaries[0], tmp_path, commands("results/testworld.root"))
    assert first.returncode == 0, first.stdout + first.stderr
    root = tmp_path / "results/testworld.root"
    cr_before, ab_before = load(root)
    second = run_macro(binaries[1], tmp_path, commands("results/lightpipe.root",
        "/bbr/thermal/emitterCenter -51 0 0 mm\n/bbr/thermal/emitterSize 1 7 7 mm"))
    assert second.returncode == 0, second.stdout + second.stderr
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


def test_output_failure_is_fatal(binaries, tmp_path):
    (tmp_path / "blocked").write_text("a file, not a directory")
    result = run_macro(binaries[0], tmp_path, commands("blocked/test.root"))
    assert result.returncode != 0, result.stdout + result.stderr
    assert "BBR022" in result.stdout + result.stderr


def test_multirun_metadata_captures_configuration_and_json_escaping(binaries, tmp_path):
    data = tmp_path / 'data"quote'
    data.mkdir()
    macro = (f"/bbr/dataDir '{data}'\n" + commands("one.root") +
             "/bbr/gun/energy_eV 0.003\n/analysis/setFileName two.root\n/run/beamOn 10\n")
    result = run_macro(binaries[0], tmp_path, macro)
    assert result.returncode == 0, result.stdout + result.stderr
    first = json.loads((tmp_path / "one.metadata.json").read_text())
    second = json.loads((tmp_path / "two.metadata.json").read_text())
    assert first["run_id"] == 0 and second["run_id"] == 1
    assert first["events"] == 20 and second["events"] == 10
    assert first["worker_count"] == 1
    assert first["configuration"] != second["configuration"]
    assert first["data"]["directory"] == str(data)
    assert first["geometry"] and first["random_engine_state"]
