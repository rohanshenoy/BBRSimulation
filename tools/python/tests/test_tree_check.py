"""The real-tree frequency check: validation/Scripts/make_tree_check.py writes a
test-world run for an HFSS data tree; validation/check_tree_run.py checks it."""
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from bbrsim import hfss

CRACK1, CRACK2, ROUND = "InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2", "RoundGap_r50um"
GRID_KEYS = {"incident_phi_deg": [0.0, 45.0, 90.0],
             "incident_theta_deg": [0.0, 45.0, 90.0, 135.0, 180.0]}
# Minimal sidecars: only the fields the generator reads (values from the legacy
# crack1 sidecar and the round-gap mock).
SLAB = {"modes": {"cutoff_ghz": 14.9896229, "polarization_filter_limit_ghz": 2997.92458},
        "exit_field": {"cross_section": {"shape": "rectangle", "y_e_half_m": 0.005, "z_e_half_m": 2.5e-05}},
        "excitation": dict(GRID_KEYS)}
DISC = {"modes": {"cutoff_ghz": 1756.984664434265, "polarization_filter_limit_ghz": 1756.984664434265},
        "exit_field": {"cross_section": {"shape": "disc", "radius_m": 5e-05}},
        "excitation": dict(GRID_KEYS)}
LEGACY = {CRACK1: ["500"], CRACK2: ["500"]}
MILESTONE = {CRACK1: ["180", "320", "569", "1010"], CRACK2: ["500"]}
GUN_KINDS = ["normal_random", "normal_E_gap", "normal_E_long", "oblique_random", "oblique_random",
             "oblique_random", "oblique_diagonal_fixed"]


def sidecar_tree(root, spec, sidecars=None):
    """A data root whose waveguides/ holds, per {id: [tokens]}, empty
    <id>_<token>GHz_Ephi={0,1} directories and a minimal sidecar each."""
    wg = Path(root) / "waveguides"
    wg.mkdir(parents=True, exist_ok=True)
    for did, tokens in spec.items():
        for tok in tokens:
            stem = f"{did}_{tok}GHz"
            (wg / f"{stem}_Ephi=0").mkdir()
            (wg / f"{stem}_Ephi=1").mkdir()
            sc = json.loads(json.dumps((sidecars or {}).get(did) or (DISC if did == ROUND else SLAB)))
            sc.update({"dataset_id": did, "frequency_label": f"{tok}GHz", "frequency_ghz": float(tok)})
            (wg / f"{stem}.dataset.json").write_text(json.dumps(sc))
    return Path(root)


def generate(repo_root, root, out, *extra):
    """(CompletedProcess, manifest or None) of make_tree_check.py root --out out."""
    r = subprocess.run([sys.executable, str(repo_root / "validation/Scripts/make_tree_check.py"),
                        str(root), "--out", str(out), *map(str, extra)], capture_output=True, text=True)
    man = json.loads((Path(out) / "tree_check.json").read_text()) if r.returncode == 0 else None
    return r, man


def kinds(man, did):
    return [r["kind"] for r in man["runs"] if r["dataset"] == did]


def macro_blocks(text):
    """The commands after /run/initialize, split at each /run/beamOn: [(commands, beamOn line)]."""
    blocks, cur = [], []
    for line in text.split("/run/initialize", 1)[1].splitlines():
        if line.startswith("/run/beamOn"):
            blocks.append((cur, line))
            cur = []
        elif line.startswith("/"):
            cur.append(line)
    return blocks


# --- the generator ----------------------------------------------------------------
def test_legacy_shaped_tree(repo_root, tmp_path):
    r, man = generate(repo_root, sidecar_tree(tmp_path / "data", LEGACY), tmp_path / "run")
    assert r.returncode == 0, r.stderr
    assert kinds(man, CRACK1) == GUN_KINDS + ["planck"] and kinds(man, CRACK2) == GUN_KINDS + ["planck"]
    assert [x["run"] for x in man["runs"]] == list(range(16))
    assert all(x["file"] == f"output/tree_r{x['run']:03d}.root" for x in man["runs"])
    assert man["not_exercised"] == [] and man["skipped"] == []
    assert all(x["expected_grid_ghz"] == 500.0 for x in man["runs"] if x["kind"] != "planck")
    assert r.stdout.strip().startswith("wrote 16 runs for 2 dataset(s) to ")


def test_milestone_shaped_tree(repo_root, tmp_path):
    r, man = generate(repo_root, sidecar_tree(tmp_path / "data", MILESTONE), tmp_path / "run")
    assert r.returncode == 0, r.stderr
    runs = man["runs"]
    assert len(runs) == 45 and sum(x["kind"].startswith("probe") for x in runs) == 8
    normal = [x for x in runs if x["dataset"] == CRACK1 and x["kind"] == "normal_random"]
    assert [x["target_ghz"] for x in normal] == pytest.approx([180 * 1.001, 320, 569, 1010 * 0.999])
    assert [x["expected_grid_ghz"] for x in normal] == [180.0, 320.0, 569.0, 1010.0]
    grid, want = [180.0, 320.0, 569.0, 1010.0], []
    for lo, hi in zip(grid[:-1], grid[1:]):
        m = math.sqrt(lo * hi)
        want += [("probe_midpoint", m * 0.99, lo), ("probe_midpoint", m * 1.01, hi)]
    want += [("probe_clamp_low", 90.0, 180.0), ("probe_clamp_high", 2020.0, 1010.0)]
    probes = [(x["kind"], x["target_ghz"], x["expected_grid_ghz"]) for x in runs if x["kind"].startswith("probe")]
    assert [(k, e) for k, _, e in probes] == [(k, e) for k, _, e in want]
    assert [t for _, t, _ in probes] == pytest.approx([t for _, t, _ in want])
    assert not any(x["kind"].startswith("probe") for x in runs if x["dataset"] == CRACK2)


def test_energies_and_expected_selection(repo_root, tmp_path):
    _, man = generate(repo_root, sidecar_tree(tmp_path / "data", MILESTONE), tmp_path / "run")
    grids = {d: [(f, s) for f, s in v["grid"]] for d, v in man["datasets"].items()}
    for x in man["runs"]:
        if x["kind"] == "planck":
            continue
        nu = hfss.photon_frequency_GHz(x["energy_eV"])
        assert nu == pytest.approx(x["target_ghz"], rel=1e-10)
        assert hfss.select_frequency(grids[x["dataset"]], nu)[0] == x["expected_grid_ghz"]


def test_directions_fold_onto_their_keys(repo_root, tmp_path):
    _, man = generate(repo_root, sidecar_tree(tmp_path / "data", LEGACY), tmp_path / "run")
    face = {CRACK1: (0.0, 0.0, 0.0), CRACK2: (0.0, 0.0, 3.0)}
    for x in man["runs"]:
        if x["kind"] == "planck":
            continue
        inc = hfss.fold_incidence(x["direction"])
        assert (inc.phi_deg, inc.theta_deg) == pytest.approx(tuple(x["incidence_key"]), abs=1e-9)
        p, k = np.array(x["position_mm"]), np.array(x["direction"])
        assert p[0] == -20.0
        assert np.allclose(p + (-p[0] / k[0]) * k, face[x["dataset"]], rtol=0, atol=1e-9)
    keys = sorted({tuple(x["incidence_key"]) for x in man["runs"]
                   if x["dataset"] == CRACK1 and x["kind"] != "planck"})
    assert keys == [(0.0, 135.0), (0.0, 180.0), (45.0, 135.0), (90.0, 135.0)]


def test_diagonal_polarization_has_equal_components(repo_root, tmp_path):
    _, man = generate(repo_root, sidecar_tree(tmp_path / "data", LEGACY), tmp_path / "run")
    diag = [x for x in man["runs"] if x["kind"] == "oblique_diagonal_fixed"]
    assert len(diag) == 2
    for x in diag:
        et, ep = hfss.polarization_components(x["polarization"], hfss.fold_incidence(x["direction"]))
        assert (et, ep) == pytest.approx((math.sqrt(0.5), math.sqrt(0.5)), abs=1e-12)
    # the vector Validation_RoundGap.mac's r08 uses at the same direction
    assert diag[0]["polarization"] == pytest.approx([0.5, 0.14644661, -0.85355339], abs=1e-8)


def test_macro_sets_every_parameter_of_every_run(repo_root, tmp_path):
    root = sidecar_tree(tmp_path / "data", MILESTONE)
    _, man = generate(repo_root, root, tmp_path / "run")
    text = (tmp_path / "run" / "tree_check.mac").read_text()
    head = text.split("/run/initialize")[0]
    assert f'/bbr/dataDir "{os.path.realpath(root)}"' in head and "/random/setSeeds 2718 28" in head
    assert "/bbr/testworld/roundGap" not in text
    blocks = macro_blocks(text)
    assert len(blocks) == len(man["runs"])
    gun = {"/bbr/gun/mode", "/bbr/gun/energy_eV", "/bbr/gun/posX", "/bbr/gun/posY", "/bbr/gun/posZ",
           "/bbr/gun/dirX", "/bbr/gun/dirY", "/bbr/gun/dirZ", "/bbr/gun/pol"}
    planck = {"/bbr/gun/mode", "/bbr/thermal/setT", "/bbr/thermal/emitterCenter", "/bbr/thermal/emitterSize"}
    for (cmds, beam), x in zip(blocks, man["runs"]):
        assert beam == f"/run/beamOn {x['events']}"
        assert f"/analysis/setFileName {x['file'][:-len('.root')]}" in cmds
        names = {c.split()[0] for c in cmds}
        if x["kind"] == "planck":
            assert planck <= names and "/bbr/gun/mode false" in cmds
        else:
            assert gun <= names and "/bbr/gun/mode true" in cmds
            e = next(c for c in cmds if c.startswith("/bbr/gun/energy_eV"))
            assert float(e.split()[1]) == x["energy_eV"]


def test_round_gap_switch_and_aim(repo_root, tmp_path):
    _, man = generate(repo_root, sidecar_tree(tmp_path / "data", {**LEGACY, ROUND: ["2000"]}), tmp_path / "run")
    text = (tmp_path / "run" / "tree_check.mac").read_text()
    assert text.index("/bbr/testworld/roundGap true") < text.index("/run/initialize")
    rg = [x for x in man["runs"] if x["dataset"] == ROUND]
    assert len(rg) == 8
    for x in rg:
        if x["kind"] != "planck":
            k = x["direction"]
            assert x["position_mm"][2] == pytest.approx(-80.0 - 20.0 / k[0] * k[2])
    p = next(x for x in rg if x["kind"] == "planck")
    assert p["emitter_center_mm"] == [-0.03, 0.0, -80.0] and p["emitter_size_mm"] == [0.02, 0.07, 0.07]


def test_unknown_dataset_not_exercised(repo_root, tmp_path):
    root = sidecar_tree(tmp_path / "data", {**LEGACY, "LightPipeSLAC": ["500"]})
    (root / "waveguides" / "LightPipeSLAC_500GHz.dataset.json").unlink()   # never read
    r, man = generate(repo_root, root, tmp_path / "run")
    assert r.returncode == 0, r.stderr
    assert man["not_exercised"] == ["LightPipeSLAC"] and man["datasets"]["LightPipeSLAC"]["aimed"] is False
    assert not any(x["dataset"] == "LightPipeSLAC" for x in man["runs"])
    assert "not exercised: LightPipeSLAC" in r.stdout


def test_grid_without_oblique_row(repo_root, tmp_path):
    flat = json.loads(json.dumps(SLAB))
    flat["excitation"]["incident_theta_deg"] = [0.0, 90.0, 180.0]
    _, man = generate(repo_root, sidecar_tree(tmp_path / "a", LEGACY, {CRACK1: flat}), tmp_path / "ra")
    assert kinds(man, CRACK1) == ["normal_random", "normal_E_gap", "normal_E_long", "planck"]
    assert man["skipped"] == [{"dataset": CRACK1, "kinds": "oblique",
                               "reason": "no declared IWaveTheta strictly between 90 and 180 deg"}]
    no45 = json.loads(json.dumps(SLAB))
    no45["excitation"]["incident_phi_deg"] = [0.0, 90.0]
    _, man = generate(repo_root, sidecar_tree(tmp_path / "b", LEGACY, {CRACK2: no45}), tmp_path / "rb")
    assert kinds(man, CRACK2).count("oblique_random") == 2
    assert "oblique_diagonal_fixed" not in kinds(man, CRACK2)
    assert man["skipped"] == [{"dataset": CRACK2, "kinds": "oblique_diagonal_fixed",
                               "reason": "IWavePhi 45 not declared (IWaveTheta 135)"}]


@pytest.mark.parametrize("case,message", [
    ("no_waveguides", "must hold waveguides/"),
    ("no_crack2", "has no dataset for InfParallelPlate_crack2"),
    ("no_sidecar", "InfParallelPlate_crack1Rohan_320GHz.dataset.json: missing"),
    ("out_inside", "refusing to write into the data tree"),
    ("legacy_dir", "legacy directory InfParallelPlate_crack2_Ephi=0"),
])
def test_refusals(repo_root, tmp_path, case, message):
    root, out = tmp_path / "data", tmp_path / "run"
    if case == "no_waveguides":
        root.mkdir()
    else:
        sidecar_tree(root, {CRACK1: ["500"]} if case == "no_crack2" else MILESTONE)
    if case == "no_sidecar":
        (root / "waveguides" / f"{CRACK1}_320GHz.dataset.json").unlink()
    if case == "out_inside":
        out = root / "run"
    if case == "legacy_dir":
        (root / "waveguides" / f"{CRACK2}_Ephi=0").mkdir()
    r, _ = generate(repo_root, root, out)
    assert r.returncode == 2 and message in r.stderr, r.stderr
    if case == "no_sidecar":
        assert "check_dataset_sidecars.py" in r.stderr
    assert not (out / "tree_check.mac").exists() and not (out / "tree_check.json").exists()
