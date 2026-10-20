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

from bbrsim import hfss, sidecar

import schema2_fixture as fx

CRACK1, CRACK2, ROUND = "InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2", "RoundGap_r50um"
GRID_KEYS = {"incident_phi_deg": [0.0, 45.0, 90.0],
             "incident_theta_deg": [0.0, 45.0, 90.0, 135.0, 180.0]}
# Minimal sidecars: only the fields the generator reads (values from the legacy
# crack1 sidecar and the round-gap mock).
SLAB = {"schema_version": "1.0",
        "modes": {"mode": "TE10", "cutoff_ghz": 14.9896229, "polarization_filter_limit_ghz": 2997.92458},
        "exit_field": {"cross_section": {"shape": "rectangle", "y_e_half_m": 0.005, "z_e_half_m": 2.5e-05}},
        "excitation": dict(GRID_KEYS)}
DISC = {"schema_version": "1.0",
        "modes": {"mode": "TE11", "cutoff_ghz": 1756.984664434265, "polarization_filter_limit_ghz": 1756.984664434265},
        "exit_field": {"cross_section": {"shape": "disc", "radius_m": 5e-05}},
        "excitation": dict(GRID_KEYS)}
LEGACY = {CRACK1: ["500"], CRACK2: ["500"]}
MILESTONE = {CRACK1: ["180", "320", "569", "1010"], CRACK2: ["500"]}
GRID35 = {"incident_phi_deg": [0.0, 22.5, 45.0, 67.5, 90.0],
          "incident_theta_deg": [90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0]}
# Minimal schema-2.0 sidecars: the fields the generator reads (crack2's section, the disc of the round gap).
SLAB2 = {"schema_version": "2.0",
         "modes": {"mode": "TE10", "cutoff_ghz": 14.9896229, "polarization_filter_limit_ghz": 1498.96},
         "cross_term": {"valid_below_ghz": 1498.96},
         "exit_field": {"cross_section": {"shape": "rectangle", "y_e_half_m": 0.005, "z_e_half_m": 5e-05}},
         "excitation": dict(GRID35)}
DISC2 = {"schema_version": "2.0",
         "modes": {"mode": "TE11", "cutoff_ghz": 1756.984664434265, "polarization_filter_limit_ghz": 1756.984664434265},
         "cross_term": {"valid_below_ghz": 1756.984664434265},
         "exit_field": {"cross_section": {"shape": "disc", "radius_m": 5e-05}},
         "excitation": dict(GRID35)}
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


def test_run_plan_35_keys(repo_root, tmp_path):
    """OR-036: the plan follows the declared 35-key grid; theta = 90 stays out of the gun runs."""
    root = sidecar_tree(tmp_path / "data", {CRACK1: ["500"], CRACK2: ["500"]},
                        {CRACK1: SLAB2, CRACK2: SLAB2})
    r, man = generate(repo_root, root, tmp_path / "run")
    assert r.returncode == 0, r.stderr
    runs = [x for x in man["runs"] if x["dataset"] == CRACK1 and x["kind"] != "planck"]
    count = lambda k: sum(x["kind"] == k for x in runs)
    assert (count("normal_random"), count("oblique_random"), count("oblique_diagonal_fixed")) == (1, 25, 5)
    assert count("normal_E_gap") == 1 and count("normal_E_long") == 1 and len(runs) == 33
    on_grid = {(p, t) for p in GRID35["incident_phi_deg"] for t in GRID35["incident_theta_deg"]}
    assert all(tuple(x["incidence_key"]) in on_grid for x in runs)
    assert not any(x["incidence_key"][1] == 90.0 for x in runs)
    assert sorted({x["incidence_key"][1] for x in runs if x["kind"] == "oblique_random"}) == \
        [105.0, 120.0, 135.0, 150.0, 165.0]
    assert sorted({x["incidence_key"][0] for x in runs if x["kind"] == "oblique_random"}) == \
        GRID35["incident_phi_deg"]
    assert [x["incidence_key"] for x in runs if x["kind"] == "oblique_diagonal_fixed"] == \
        [[45.0, t] for t in (105.0, 120.0, 135.0, 150.0, 165.0)]


def test_manifest_v2_fields(repo_root, tmp_path):
    """The manifest is /2 and records, per dataset, the schema major, the regime cutoffs the log is
    checked against and the cross-term limit (null under schema 1)."""
    root = sidecar_tree(tmp_path / "data", {CRACK1: ["500"], CRACK2: ["500"], ROUND: ["2000"]},
                        {CRACK1: SLAB2, CRACK2: SLAB2, ROUND: DISC2})
    r, man = generate(repo_root, root, tmp_path / "run")
    assert r.returncode == 0, r.stderr
    assert man["schema"] == "bbrsim-tree-check/2"
    c2 = man["datasets"][CRACK2]
    assert c2["schema_major"] == 2 and c2["rho_valid_below_ghz"] == pytest.approx(1498.96)
    modes = [m for m, _ in c2["regime_cutoffs"]]
    assert modes[0] == "TE10" and "TE01" in modes and len(modes) == 14      # TE10 + the gap family, 0.1 mm
    assert [g for _, g in c2["regime_cutoffs"]] == sorted(g for _, g in c2["regime_cutoffs"])
    assert dict(c2["regime_cutoffs"])["TE01"] == pytest.approx(1498.962, rel=1e-6)
    rg = man["datasets"][ROUND]
    assert [m for m, _ in rg["regime_cutoffs"]] == ["TE11", "TM01", "TE21"]
    assert rg["regime_cutoffs"][0][1] == pytest.approx(1756.9847, rel=1e-6) and rg["schema_major"] == 2
    r, man = generate(repo_root, sidecar_tree(tmp_path / "data1", LEGACY), tmp_path / "run1")
    assert r.returncode == 0, r.stderr
    d1 = man["datasets"][CRACK1]
    assert man["schema"] == "bbrsim-tree-check/2" and d1["schema_major"] == 1 and d1["rho_valid_below_ghz"] is None
    assert d1["regime_cutoffs"] == [["TE10", pytest.approx(14.9896229)]]    # schema 1: the lowest mode only


# --- the validator: pure functions -------------------------------------------------
@pytest.fixture(scope="module")
def checker(repo_root):
    spec = importlib.util.spec_from_file_location("check_tree_run", repo_root / "validation/check_tree_run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_expected_warnings(checker):
    grids = {"a": [(180.0, "a_180GHz"), (320.0, "a_320GHz")], "b": [(500.0, "b_500GHz")]}
    cutoffs = {"a": [("TE10", 14.9896229)], "b": [("TE10", 14.9896229)]}
    none = {"a": None, "b": None}
    want = checker.expected_warnings([("a", 200.0, 180.0), ("b", 1e4, 500.0)], grids, cutoffs, none)
    assert set(want.values()) == {0}            # inside the grid; a one-point grid never clamps
    want = checker.expected_warnings([("a", 100.0, 180.0), ("a", 2000.0, 320.0), ("a", 12.0, 180.0),
                                      ("b", 12.0, 500.0)], grids, cutoffs, none)
    assert want[("BBR008", "a", "low")] == 1 and want[("BBR008", "a", "high")] == 1
    assert want[("BBR026", "a", "TE10", "below")] == 1 and want[("BBR026", "b", "TE10", "below")] == 1
    assert want[("BBR008", "b", "low")] == 0 and want[("BBR026", "a", "TE10", "above")] == 0
    # 'above': a photon at or above the cutoff served from a grid point below it
    fc = 1756.984664434265
    disc, cut = {"r": [(1000.0, "r_1000GHz"), (2000.0, "r_2000GHz")]}, {"r": [("TE11", fc)]}
    want = checker.expected_warnings([("r", 1400.0, 1000.0), ("r", fc, 1000.0)], disc, cut, {"r": None})
    assert want[("BBR026", "r", "TE11", "above")] == 1 and want[("BBR026", "r", "TE11", "below")] == 0


def test_expected_warnings_per_mode(checker):
    d = CRACK2
    grids = {d: [(500.0, f"{d}_500GHz")]}
    cutoffs = {d: [("TE10", 14.9896229), ("TE01", 1498.96)]}
    want = checker.expected_warnings([(d, 10.0, 500.0), (d, 1600.0, 500.0)], grids, cutoffs, {d: None})
    assert {k: v for k, v in want.items() if k[0] == "BBR026"} == {
        ("BBR026", d, "TE10", "below"): 1, ("BBR026", d, "TE10", "above"): 0,
        ("BBR026", d, "TE01", "below"): 0, ("BBR026", d, "TE01", "above"): 1}
    # a photon between the two cutoffs, served from the 500 GHz table, crosses neither
    want = checker.expected_warnings([(d, 100.0, 500.0)], grids, cutoffs, {d: None})
    assert [k for k, v in want.items() if v] == []


def test_expected_bbr028(checker):
    fc = 1756.984664434265
    grids = {"r": [(2000.0, "r_2000GHz")], "s": [(2000.0, "s_2000GHz")]}
    cutoffs = {"r": [("TE11", fc)], "s": [("TE11", fc)]}
    want = checker.expected_warnings([("r", 2000.0, 2000.0), ("r", 1900.0, 2000.0), ("s", 2000.0, 2000.0)],
                                     grids, cutoffs, {"r": fc, "s": None})
    assert want[("BBR028", "r")] == 1             # once per dataset, however many entries
    assert want[("BBR028", "s")] == 0             # a schema-1 dataset never warns
    below = checker.expected_warnings([("r", 1000.0, 1000.0)], {"r": [(1000.0, "r_1000GHz")]},
                                      {"r": [("TE11", fc)]}, {"r": fc})
    assert below[("BBR028", "r")] == 0            # a table below the limit does not warn
    # the grid point decides, not the photon: a photon below the limit served from a table above it
    edge = checker.expected_warnings([("r", 1700.0, 2000.0)], {"r": [(2000.0, "r_2000GHz")]},
                                     {"r": [("TE11", fc)]}, {"r": fc})
    assert edge[("BBR028", "r")] == 1


LOG = """[BBR] HFSS dataset InfParallelPlate_crack1Rohan: 1 frequency grid point(s): 500 GHz (from /t/waveguides)
[BBR] crack InfParallelPlate_crack1Rohan: HFSS (p, l, g) = (1, 10, 0.05) mm, Geant4 (x, y, z) = (1, 10, 0.05) mm; 1 sidecar(s) fit
G4WT1 > 
-------- WWWW ------- G4Exception-START -------- WWWW -------
*** G4Exception : BBR026
      issued by : BBRCrackLibrary::Lookup
BBR026 dataset=InfParallelPlate_crack1Rohan direction=below mode=TE10 nu_GHz=11.4 grid_GHz=500 cutoff_GHz=14.9896 : the photon and the HFSS dataset serving it lie on opposite sides of a regime cutoff; the table is used unchanged. Reported once per dataset, mode and direction.
*** This is just a warning message. ***
-------- WWWW -------- G4Exception-END --------- WWWW -------
G4WT3 > [BBR] HFSS InfParallelPlate_crack1Rohan_500GHz key (0, 180): max transmittance 1.05451 > 1 (port-normalization artefact), normalized to 1
"""


def test_parse_log(checker):
    lg = checker.parse_log(LOG)
    assert lg["warnings"] == {("BBR026", CRACK1, "TE10", "below"): 1}
    assert lg["banners"] == 1 and lg["geomnav"] == 0 and lg["other_codes"] == []
    assert len(lg["crack_lines"]) == 1 and len(lg["dataset_lines"]) == 1
    bad = checker.parse_log(LOG + "G4Exception-START\n*** G4Exception : BBR010\nGeomNav1002: stuck track\n")
    assert len(bad["other_codes"]) == 1 and bad["banners"] == 2 and bad["geomnav"] == 1


LOG2 = """[BBR] crack RoundGap_r50um: HFSS (p, l, g) = (0.4, 0.1, 0.1) mm, Geant4 (x, y, z) = (0.4, 0.1, 0.1) mm; 1 sidecar(s) fit
-------- WWWW ------- G4Exception-START -------- WWWW -------
*** G4Exception : BBR026
      issued by : BBRCrackLibrary::Lookup
BBR026 dataset=InfParallelPlate_crack2 direction=above mode=TE01 nu_GHz=1600 grid_GHz=500 cutoff_GHz=1498.96 : the table is used unchanged.
*** This is just a warning message. ***
-------- WWWW -------- G4Exception-END --------- WWWW -------
-------- WWWW ------- G4Exception-START -------- WWWW -------
*** G4Exception : BBR028
      issued by : BBRCrackLibrary::Lookup
BBR028 dataset=RoundGap_r50um grid_GHz=2000 valid_below_GHz=1756.98 : the cross term is not valid above the polarization-filter limit.
*** This is just a warning message. ***
-------- WWWW -------- G4Exception-END --------- WWWW -------
"""


def test_parse_log_bbr026_mode_and_bbr028(checker):
    lg = checker.parse_log(LOG2)
    assert lg["warnings"] == {("BBR026", CRACK2, "TE01", "above"): 1, ("BBR028", ROUND): 1}
    assert lg["banners"] == 2 and lg["other_codes"] == []        # OR-024: BBR028 is whitelisted
    lg = checker.parse_log(LOG2 + "BBR029 something else\n")
    assert len(lg["other_codes"]) == 1


def test_planck_prediction_clamped_s2(checker):
    """CDX-016: under 2.x the Planck run's per-photon T is the clamped mixture mean, not (T0 + T1)/2."""
    far = dict(theta_deg=np.array([90.0]), phi_deg=np.array([0.0]), F0_theta=np.array([1 + 0j]),
               F0_phi=np.array([0j]), F1_theta=np.array([0j]), F1_phi=np.array([1 + 0j]))
    a = hfss.AngleDataset(key=(0.0, 180.0), T0=1.6, T1=0.2, schema=2, **far)
    t = checker.random_polarization_T(a)
    assert t == pytest.approx(hfss.random_polarization_mixture(a)[0], rel=1e-14)
    assert t < 0.5 * (1.6 + 0.2) - 0.01 and 0.0 < t < 1.0
    # the closed form: mean over psi of min(1, 1.6 cos^2 + 0.2 sin^2)
    psi = (np.arange(360) + 0.5) / 360 * 2 * np.pi
    assert t == pytest.approx(np.minimum(1.0, 1.6 * np.cos(psi) ** 2 + 0.2 * np.sin(psi) ** 2).mean(), rel=1e-12)
    # schema 1 keeps (T0 + T1)/2 exactly
    b = hfss.AngleDataset(key=(0.0, 180.0), T0=0.6, T1=0.3, schema=1, **far)
    assert checker.random_polarization_T(b) == pytest.approx(0.45, rel=1e-14)


def test_manifest_v1_refused(checker, tmp_path):
    root = sidecar_tree(tmp_path / "data", LEGACY)
    run = tmp_path / "run"
    run.mkdir()
    (run / "tree_check.json").write_text(json.dumps({"schema": "bbrsim-tree-check/1", "data_root": str(root),
                                                     "datasets": {}, "runs": []}))
    rows = checker.validate(str(run)).rows
    assert rows[0][0] == "FAIL" and rows[0][1] == "manifest" and "bbrsim-tree-check/2" in rows[0][2]
    assert len(rows) == 1


def tree_2x(tmp_path):
    """A one-dataset 2.x tree ('crack', 500 GHz, the nine-key fixture) with a complete sidecar, and its
    manifest in a run directory."""
    base = tmp_path / "data" / "waveguides"
    base.mkdir(parents=True)
    stem = fx.write_nine_key_grid2(str(base))
    sc = fx.sidecar2()
    sc["files"] = sidecar.file_entries(str(base), stem)
    sidecar.write(str(base), stem, sc)
    sha = hashlib_sha(base / f"{stem}.dataset.json")
    man = {"schema": "bbrsim-tree-check/2", "data_root": str(tmp_path / "data"),
           "datasets": {"crack": {"aimed": True, "grid": [[500.0, stem]], "schema_major": 2,
                                  "cutoff_ghz": 16.65, "regime_cutoffs": [["TE10", 16.65]],
                                  "rho_valid_below_ghz": 1000.0, "sidecar_sha256": {stem: sha}}},
           "not_exercised": [], "skipped": [], "runs": []}
    run = tmp_path / "run"
    run.mkdir()
    (run / "tree_check.json").write_text(json.dumps(man))
    return run, base, stem


def hashlib_sha(path):
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_tree_check_runs_check_full_on_2x_sidecars(checker, tmp_path):
    """The tables of a 2.x tree are loaded by hfss.load_dataset, which does not repeat the 2.x CSV row
    checks; the sidecar row runs sidecar.check_full (checksums, InterceptFactor, T sources, dark rows)."""
    run, base, stem = tree_2x(tmp_path)
    row = next(r for r in checker.validate(str(run)).rows if r[1] == "sidecars 2.x")
    assert row[0] == "PASS", row
    # an InterceptFactor edited in one row, checksums refreshed: only check_full sees it
    path = base / f"{stem}_Ephi=0" / "waveguide.csv"
    lines = path.read_text().splitlines()
    v = lines[1].split(",")
    v[-1] = repr(float(v[-1]) * 1.01)
    lines[1] = ",".join(v)
    path.write_text("\n".join(lines) + "\n")
    sc = sidecar.load(str(base), stem)
    sc["files"] = sidecar.file_entries(str(base), stem)
    sidecar.write(str(base), stem, sc)
    row = next(r for r in checker.validate(str(run)).rows if r[1] == "sidecars 2.x")
    assert row[0] == "FAIL" and "BBR012" in row[2] and stem in row[2], row


def test_no_sidecar_row_under_schema_1(checker, tmp_path):
    root = sidecar_tree(tmp_path / "data", LEGACY)
    man = {"schema": "bbrsim-tree-check/2", "data_root": str(root),
           "datasets": {d: {"aimed": True, "grid": [[500.0, f"{d}_500GHz"]], "schema_major": 1,
                            "cutoff_ghz": 14.98, "regime_cutoffs": [["TE10", 14.98]],
                            "rho_valid_below_ghz": None} for d in (CRACK1, CRACK2)},
           "not_exercised": [], "skipped": [], "runs": []}
    run = tmp_path / "run"
    run.mkdir()
    (run / "tree_check.json").write_text(json.dumps(man))
    assert not any(n.startswith("sidecars") for _, n, _ in checker.validate(str(run)).rows)


def test_z_max(checker):
    assert [round(checker.z_max(m), 4) for m in (0, 1, 15, 16, 300, 5000)] == \
        [4.0, 4.0, 4.0, 4.0032, 4.6491, 5.1993]


def test_report_family_wise(checker):
    for off, passes in ((4.5, True), (10.0, False)):          # the limit is 4.65 SE at M = 300
        rep = checker.Report()
        for i in range(299):
            rep.stat(f"row {i}", 100.0, 100.0, 1.0)
        rep.stat("off", 100.0 + off, 100.0, 1.0)
        assert rep.finish() is passes
    rep = checker.Report()
    rep.stat("one short at T = 1", 299.0, 300.0, 0.0, slack=1.0)   # one count of slack when sigma is 0
    assert rep.finish() is True
    rep = checker.Report()
    rep.stat("two short at T = 1", 298.0, 300.0, 0.0, slack=1.0)
    assert rep.finish() is False


def test_validator_skip_lines(checker, tmp_path):
    root = sidecar_tree(tmp_path / "data", {**LEGACY, "LightPipeSLAC": ["500"]})
    man = {"schema": "bbrsim-tree-check/2", "data_root": str(root),
           "datasets": {d: {"aimed": d != "LightPipeSLAC", "grid": [[500.0, f"{d}_500GHz"]],
                            "cutoff_ghz": 14.9896229, "schema_major": 1,
                            "regime_cutoffs": [["TE10", 14.9896229]], "rho_valid_below_ghz": None,
                            "sidecar_sha256": {f"{d}_500GHz": checker.sidecar_sha256(str(root / "waveguides"),
                                                                                    f"{d}_500GHz")}}
                     for d in (CRACK1, CRACK2, "LightPipeSLAC")},
           "not_exercised": ["LightPipeSLAC"],
           "skipped": [{"dataset": CRACK1, "kinds": "oblique", "reason": "no row"}], "runs": []}
    run = tmp_path / "run"
    run.mkdir()
    (run / "tree_check.json").write_text(json.dumps(man))
    rows = checker.validate(str(run)).rows
    assert ("SKIP", "dataset LightPipeSLAC", "the test world cannot place it; not exercised") in rows
    assert ("SKIP", f"dataset {CRACK1} oblique", "no row") in rows
    assert any(s == "PASS" and n == "tree unchanged" for s, n, _ in rows)
    assert any(s == "FAIL" and n == "selection" for s, n, _ in rows)      # no run, no entry: fails closed


def test_tree_changed_after_generation(repo_root, checker, tmp_path):
    root = sidecar_tree(tmp_path / "data", LEGACY)
    r, _ = generate(repo_root, root, tmp_path / "run")
    assert r.returncode == 0, r.stderr
    (root / "waveguides" / f"{CRACK1}_320GHz_Ephi=0").mkdir()        # a frequency added after generation
    rows = checker.validate(str(tmp_path / "run")).rows
    assert next(s for s, n, _ in rows if n == "tree unchanged") == "FAIL"


def test_tree_rebuilt_with_same_labels(repo_root, checker, tmp_path):
    # A tree rebuilt in place keeps its labels but not its sidecars, which carry the CSV checksums.
    root = sidecar_tree(tmp_path / "data", LEGACY)
    r, _ = generate(repo_root, root, tmp_path / "run")
    assert r.returncode == 0, r.stderr
    sc = root / "waveguides" / f"{CRACK1}_500GHz.dataset.json"
    sc.write_text(sc.read_text().replace("}", ', "rebuilt": true}', 1))
    rows = checker.validate(str(tmp_path / "run")).rows
    assert next(s for s, n, _ in rows if n == "tree unchanged") == "FAIL"


# --- the exit-direction row (CDX-027) ---------------------------------------------------
D_OUT = (0.0010000005, -0.7071067, 0.7071067)


def crossing(nb, vol_pre, mat_pre, vol_post, mat_post, d_pre, d_post, status="BBRNone", x=2.0, event=0, track=1):
    return {"event_id": event, "track_id": track, "n_boundary": nb, "status": status,
            "vol_pre": vol_pre, "mat_pre": mat_pre, "vol_post": vol_post, "mat_post": mat_post,
            "px_pre": d_pre[0], "py_pre": d_pre[1], "pz_pre": d_pre[2],
            "px_post": d_post[0], "py_post": d_post[1], "pz_post": d_post[2], "x": x, "y": 0.0, "z": 0.0}


def entry(nb=1, d=D_OUT, **kw):
    return crossing(nb, "World", "G4_Galactic", "crack1", "vacuum_wg", (1.0, 0.0, 0.0), d,
                    status="BBRDiffractionTransmit", **kw)


def leave(nb=2, d_pre=D_OUT, d_post=D_OUT, **kw):
    return crossing(nb, "crack1", "vacuum_wg", "World", "G4_Galactic", d_pre, d_post, **kw)


def frame(*rows):
    import pandas as pd
    return pd.DataFrame(list(rows))


def exit_row(checker, *rows):
    rep = checker.Report()
    checker.check_exit_direction(rep, "run 000", frame(*rows))
    (status, name, detail), = rep.rows
    return status, name, detail


def test_exit_direction_clean_pair_passes(checker):
    status, name, detail = exit_row(checker, entry(), leave())
    assert (status, name) == ("PASS", "run 000 exit direction"), detail
    assert detail.startswith("1 transmitted photons")


def test_exit_direction_no_transmissions_passes_with_a_note(checker):
    status, _, detail = exit_row(checker, crossing(1, "World", "G4_Galactic", "CuPlate", "Cu_RRR100_T4K",
                                                   (1, 0, 0), (-1, 0, 0), status="BBRReflect"))
    assert status == "PASS" and "no transmissions" in detail


def test_exit_direction_wall_hit_fails_even_with_the_pre_direction_intact(checker):
    # At a wall hit the photon reaches the wall along the sampled direction, so px_pre still equals it;
    # only the volume and the reflected post direction give the event away.
    refl = (-D_OUT[0], D_OUT[1], D_OUT[2])
    status, _, detail = exit_row(checker, entry(), crossing(2, "crack1", "vacuum_wg", "CuPlate1", "Cu_RRR100_T4K",
                                                           D_OUT, refl))
    assert status == "FAIL" and "1 wall hit" in detail and "event 0 track 1 n_boundary 1" in detail
    assert "direction mismatch" not in detail


def test_exit_direction_other_volume_fails(checker):
    status, _, detail = exit_row(checker, entry(), crossing(2, "crack1", "vacuum_wg", "Detector", "G4_Si",
                                                           D_OUT, D_OUT))
    assert status == "FAIL" and "1 other volume" in detail and "wall hit" not in detail


def test_exit_direction_direction_tolerance(checker):
    off = lambda e: (D_OUT[0] + e, D_OUT[1], D_OUT[2])
    assert exit_row(checker, entry(), leave(d_post=off(1e-12), d_pre=off(1e-12)))[0] == "PASS"
    for kw in ({"d_post": off(1e-6)}, {"d_pre": off(1e-6)}):
        status, _, detail = exit_row(checker, entry(), leave(**kw))
        assert status == "FAIL" and "1 direction mismatch" in detail


def test_exit_direction_missing_following_row_fails(checker):
    status, _, detail = exit_row(checker, entry())
    assert status == "FAIL" and "1 no following row" in detail


def test_exit_direction_follows_by_track_and_counts_kinds(checker):
    # Two events interleaved, one clean, one wall hit, one more with no row: three kinds counted apart.
    rows = [entry(event=0), entry(event=1), entry(event=2), leave(event=0),
            crossing(2, "crack1", "vacuum_wg", "CuPlate1", "Cu_RRR100_T4K", D_OUT, (-D_OUT[0],) + D_OUT[1:], event=1)]
    status, _, detail = exit_row(checker, *rows)
    assert status == "FAIL" and "1 wall hit" in detail and "1 no following row" in detail
    assert "event 1 " in detail and "event 2 " in detail and "of 3 transmitted" in detail


def test_exit_direction_coplanar_plate_pair_is_a_clean_exit(checker):
    status, _, detail = exit_row(
        checker, entry(),
        crossing(2, "crack1", "vacuum_wg", "CuPlate1", "Cu_RRR100_T4K", D_OUT, D_OUT, x=2.0),
        crossing(3, "CuPlate1", "Cu_RRR100_T4K", "World", "G4_Galactic", D_OUT, D_OUT, x=2.0))
    assert status == "PASS" and "1 through coplanar" in detail
    # a real passage through the plate (a longer step) is a wall hit
    status, _, detail = exit_row(
        checker, entry(),
        crossing(2, "crack1", "vacuum_wg", "CuPlate1", "Cu_RRR100_T4K", D_OUT, D_OUT, x=2.0),
        crossing(3, "CuPlate1", "Cu_RRR100_T4K", "World", "G4_Galactic", D_OUT, D_OUT, x=2.5))
    assert status == "FAIL" and "1 wall hit" in detail


# --- the validator end to end: a small real run of the test world --------------------
@pytest.fixture(scope="module")
def e2e(repo_root, tmp_path_factory):
    """A small generated run of the test world on BBRSIMDATA's tree, with its log."""
    binary, data = os.environ.get("BBR_TESTWORLD_BINARY"), os.environ.get("BBRSIMDATA")
    if not binary or not data:
        pytest.skip("needs BBR_TESTWORLD_BINARY and BBRSIMDATA")
    run = tmp_path_factory.mktemp("tree_e2e")
    r, man = generate(repo_root, data, run, "--events", 600, "--planck-events", 3000)
    assert r.returncode == 0, r.stderr
    env = dict(os.environ, G4FORCENUMBEROFTHREADS="1")
    env.pop("DYLD_LIBRARY_PATH", None)
    env.pop("LD_LIBRARY_PATH", None)
    with open(run / "run.log", "w") as log:
        p = subprocess.run([binary, "tree_check.mac"], cwd=run, env=env, stdout=log,
                           stderr=subprocess.STDOUT, timeout=900)
    assert p.returncode == 0, (run / "run.log").read_text()[-3000:]
    return run, man


@pytest.fixture(scope="module")
def e2e_tables(checker, e2e):
    """One table cache for every end-to-end test (each 500 GHz table loads once)."""
    return checker.Tables(os.path.join(e2e[1]["data_root"], "waveguides"))


def rows_of(checker, run, tables):
    rep = checker.validate(str(run), str(run / "run.log"), tables=tables)
    return rep.finish(), rep.rows


def failing(rows):
    return {n: d for s, n, d in rows if s == "FAIL"}


def copy_run(run, dst):
    shutil.copytree(run, dst)
    return dst


def rewrite_crossings(path, change, ntuple="crossings"):
    """Rewrite the ROOT file at path with change(arrays of ntuple) applied in place."""
    import uproot
    with uproot.open(path) as f:
        trees = {name: {k: np.array(v) for k, v in f[name].arrays(library="np").items()}
                 for name in ("crossings", "abspoints")}                 # writable copies
    change(trees[ntuple])
    with uproot.recreate(path) as f:
        for name, arrs in trees.items():
            tree = f.mktree(name, {k: v.dtype for k, v in arrs.items()})
            if len(next(iter(arrs.values()))):
                tree.extend(arrs)


def test_e2e_passes(checker, e2e, e2e_tables):
    run, _ = e2e
    ok, rows = rows_of(checker, run, e2e_tables)
    assert ok, failing(rows)
    names = {n for _, n, _ in rows}
    assert {"tree unchanged", "sentinel", "selection", f"coverage {CRACK1}", f"coverage {CRACK2}",
            f"log startup {CRACK1}", f"log grid {CRACK2}", "log hygiene",
            f"log BBR008 {CRACK1} low", f"log BBR026 {CRACK2} TE10 below"} <= names
    passing = {n for s, n, _ in rows if s == "PASS"}
    assert any(n.endswith(" event_id") for n in passing)
    assert any(n.endswith(" exit direction") for n in passing)


def test_e2e_cli(repo_root, e2e):
    run, _ = e2e
    r = subprocess.run([sys.executable, str(repo_root / "validation/check_tree_run.py"), str(run),
                        "--log", str(run / "run.log")], capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.rstrip().endswith("RESULT: PASS"), r.stdout[-3000:]
    assert not r.stderr


def test_e2e_wrong_frequency(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")

    def change(cr):
        i = int(np.flatnonzero(cr["hfss_freq_GHz"] > 0)[0])
        cr["hfss_freq_GHz"][i] = 501.0

    rewrite_crossings(run / "output/tree_r000.root", change)
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and "selection" in failing(rows)


def test_e2e_event_id_out_of_range(checker, e2e, e2e_tables, tmp_path):
    run, man = e2e
    run = copy_run(run, tmp_path / "run")
    events = man["runs"][0]["events"]

    def change(cr):
        cr["event_id"][0] = events          # the first valid id is 0, the last events - 1

    rewrite_crossings(run / "output/tree_r000.root", change)
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and "run 000 event_id" in failing(rows)
    detail = failing(rows)["run 000 event_id"]
    assert "tree_r000.root" in detail and "crossings" in detail and str(events) in detail
    assert not any(n.endswith(" event_id") for n in failing(rows) if n != "run 000 event_id")


def test_e2e_event_id_negative_in_abspoints(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")

    def change(ab):
        assert len(ab["event_id"]), "run 0 has no abspoints rows to corrupt"
        ab["event_id"][-1] = -1

    rewrite_crossings(run / "output/tree_r000.root", change, ntuple="abspoints")
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and "abspoints" in failing(rows).get("run 000 event_id", "")


def test_e2e_metadata_other_tree(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")
    path = run / "output/tree_r003.metadata.json"
    meta = json.loads(path.read_text())
    meta["data"]["directory"] = "/somewhere/else"
    path.write_text(json.dumps(meta))
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and any(n.startswith("run 003 ") and n.endswith(" metadata") for n in failing(rows))


def test_e2e_extra_warning_in_log(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")
    with open(run / "run.log", "a") as fh:
        fh.write("-------- WWWW ------- G4Exception-START -------- WWWW -------\n*** G4Exception : BBR008\n"
                 f"BBR008 dataset={CRACK1} side=high nu_GHz=600 edge_GHz=500 : photon frequency lies outside\n")
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and f"log BBR008 {CRACK1} high" in failing(rows)


def test_e2e_missing_output(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")
    (run / "output/tree_r005.root").unlink()
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and any(n.startswith("run 005 ") and n.endswith(" output") for n in failing(rows))


def test_e2e_empty_planck_run(checker, e2e, e2e_tables, tmp_path):
    # Without the log: a Planck run with no crack entries must fail on its own rows.
    run, man = e2e
    run = copy_run(run, tmp_path / "run")
    r = next(x for x in man["runs"] if x["kind"] == "planck" and x["dataset"] == CRACK1)
    rewrite_crossings(run / r["file"], lambda c: c.update({k: v[:0] for k, v in c.items()}))
    rep = checker.validate(str(run), tables=e2e_tables)
    assert not rep.finish()
    assert f"{checker.run_tag(r)} entries" in failing(rep.rows)


def test_e2e_metadata_events_differ(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")
    meta = run / "output/tree_r000.metadata.json"
    m = json.loads(meta.read_text())
    m["events"] = int(m["events"]) - 1
    meta.write_text(json.dumps(m))
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and any(n.startswith("run 000 ") and n.endswith(" metadata") for n in failing(rows))


# --- the validator: statistics and the physics row ------------------------------------
FF_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag\n"
WG_HDR = ("Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,"
          "Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag\n")


def write_table(wg, stem, T, label="500GHz"):
    """A one-key (0, 180) dataset with T = (T0, T1), a two-row far field and one exit point."""
    for e, t in enumerate(T):
        d = Path(wg) / f"{stem}_Ephi={e}"
        d.mkdir(parents=True)
        (d / "far_field.csv").write_text(FF_HDR + "".join(
            f"{label},{e},0,180,0,{th},0.1,0.0,0.2,0.05\n" for th in (60, 120)))
        (d / "waveguide.csv").write_text(WG_HDR + f"{label},{e},0,180,{t},1.0,0,0,0,0,1,0,0,0,0\n")


def test_planck_fractions(checker):
    p = checker.planck_fractions([180.0, 320.0, 569.0, 1010.0], 20.0, (4.14e-5, 8.27e-2))
    assert p.sum() == pytest.approx(1.0, abs=1e-12) and (p > 0).all()
    from bbrsim import physics
    e = np.linspace(4.14e-5, 8.27e-2, 400001)
    pdf = physics.planck_photon_number_pdf(e, 20.0)
    below = np.trapezoid(np.where(hfss.photon_frequency_GHz(e) < math.sqrt(180.0 * 320.0), pdf, 0.0), e)
    assert p[0] == pytest.approx(below / np.trapezoid(pdf, e), rel=1e-3)


def test_physics_row_swapped_ephi(checker, tmp_path):
    wg = tmp_path / "waveguides"
    write_table(wg, "good_500GHz", (0.99, 1e-10))
    write_table(wg, "swap_500GHz", (1e-10, 0.99))          # Ephi directories swapped
    man = {"datasets": {d: {"aimed": True, "cross_section": "rectangle",
                            "polarization_filter_limit_ghz": 2997.92458} for d in ("good", "swap")}}
    grids = {d: [(500.0, f"{d}_500GHz")] for d in ("good", "swap")}
    rep = checker.Report()
    checker.physics_rows(rep, man, checker.Tables(str(wg)), grids)
    assert {n: s for s, n, _ in rep.rows} == {"physics good 500 GHz": "PASS", "physics swap 500 GHz": "FAIL"}
    man["datasets"]["good"]["cross_section"] = "disc"
    man["datasets"]["swap"]["polarization_filter_limit_ghz"] = 510.0     # 500 GHz is above 0.95 x 510
    rep = checker.Report()
    checker.physics_rows(rep, man, checker.Tables(str(wg)), grids)
    assert [s for s, _, _ in rep.rows] == ["SKIP", "SKIP"]


def test_e2e_statistical_rows(checker, e2e, e2e_tables):
    run, man = e2e
    ok, rows = rows_of(checker, run, e2e_tables)
    assert ok, failing(rows)
    names = [n for _, n, _ in rows]
    assert sum(n.endswith(" transmittance") for n in names) == sum(x["kind"] != "planck" for x in man["runs"])
    assert any(n.endswith(" <py^2>") for n in names)
    assert any(" transmittance InfParallelPlate_crack2 500 GHz" in n for n in names)
    assert ("PASS", f"physics {CRACK1} 500 GHz") in [(s, n) for s, n, _ in rows]


def test_e2e_transmittance_off(checker, e2e, e2e_tables, tmp_path):
    run = copy_run(e2e[0], tmp_path / "run")
    legend = json.loads((run / "output/tree_r001.metadata.json").read_text())["legend"]["status"]
    code = {name: int(c) for c, name in legend.items()}

    def change(cr):        # run 1, E across the gap at normal incidence: every transmission now reflects
        cr["status_code"][cr["status_code"] == code["BBRDiffractionTransmit"]] = code["BBRDiffractionReflect"]

    rewrite_crossings(run / "output/tree_r001.root", change)
    ok, rows = rows_of(checker, run, e2e_tables)
    assert not ok and any(n.startswith("run 001 ") and n.endswith(" transmittance") for n in failing(rows))
