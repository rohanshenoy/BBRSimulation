"""Validator negative tests and the mock HFSS generator, run as subprocesses.

check_invariants.py runs on synthetic outputs from
validation/Scripts/tests/make_bad_output.py: a clean file passes, and each seeded
violation fails, including a photon starting in copper when the legend lacks
the copper code (a NaN decode must not pass as "not a metal").
"""
import os
import subprocess
import sys

import numpy as np
import pytest

from bbrsim import hfss
from bbrsim import sidecar


def run(script, *args, cwd):
    return subprocess.run([sys.executable, str(script), *map(str, args)],
                          capture_output=True, text=True, cwd=cwd)


@pytest.mark.parametrize("switches,passes,signature", [
    ([], True, "RESULT: PASS"),
    (["--metal-start"], False, "rows STARTING in a metal  : 1"),
    (["--metal-start", "--drop-legend-code", "material:0"], False, "legend lacks code(s) [0] for mat_pre"),
    (["--unknown-term"], False, "1. 'unknown' labels                     : 1"),
    (["--empty"], False, "empty output cannot pass"),
], ids=["clean", "metal-start", "metal-start-legend-gap", "unknown-term", "empty"])
def test_check_invariants(repo_root, tmp_path, switches, passes, signature):
    mk = run(repo_root / "validation/Scripts/tests/make_bad_output.py", tmp_path / "out", *switches,
             cwd=tmp_path)
    assert mk.returncode == 0, mk.stdout + mk.stderr
    r = run(repo_root / "validation/check_invariants.py", tmp_path / "out" / "bbr.root", cwd=tmp_path)
    assert (r.returncode == 0) == passes, r.stdout + r.stderr
    assert ("RESULT: PASS" if passes else "RESULT: FAIL") in r.stdout
    assert signature in r.stdout


def test_check_invariants_empty_allowed(repo_root, tmp_path):
    # The FAIL above is the emptiness: with the flag the same file passes.
    run(repo_root / "validation/Scripts/tests/make_bad_output.py", tmp_path / "out", "--empty", cwd=tmp_path)
    r = run(repo_root / "validation/check_invariants.py", tmp_path / "out" / "bbr.root",
            "--allow-no-crossings", cwd=tmp_path)
    assert r.returncode == 0 and "RESULT: PASS" in r.stdout, r.stdout


FF_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag\n"
WG_HDR = ("Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,"
          "Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag\n")


def write_source(waveguides, dataset_id):
    """One key (0, 180), far field at Theta 60 and 120 deg, T0 = 0.5, T1 = 0.25."""
    for e, out in ((0, 0.5), (1, 0.25)):
        d = waveguides / f"{dataset_id}_500GHz_Ephi={e}"
        d.mkdir(parents=True)
        (d / "far_field.csv").write_text(FF_HDR + "".join(
            f"500GHz,{e},0,180,0,{t},0.1,0.0,0.2,0.05\n" for t in (60, 120)))
        (d / "waveguide.csv").write_text(WG_HDR + f"500GHz,{e},0,180,{out},1.0,0,0,0,0,1,0,0,0,0\n")
    sc = sidecar.build_from_csvs(
        str(waveguides), dataset_id, "500GHz",
        section={"shape": "rectangle", "y_e_half_m": 2e-3, "z_e_half_m": 2.5e-5},
        extent_mm={"p": 1.0, "l": 4.0, "g": 0.05}, provenance={"producer": "test"},
        exit_origin_mm=[0.0, 0.0, 1.0], plane_wave_origin_mm=[0.0, 0.0, 0.0],
        bounding_box_mm=[-0.025, -2.0, 0.0, 0.025, 2.0, 1.0])
    sidecar.write(str(waveguides), f"{dataset_id}_500GHz", sc)


def test_mock_generator(repo_root, tmp_path):
    gen = repo_root / "validation/Scripts/make_mock_hfss_frequencies.py"
    src = tmp_path / "tree" / "data" / "waveguides"
    write_source(src, "c")
    # refuses to write inside data/: rc 2, nothing written
    inside = tmp_path / "tree" / "data" / "mock"
    r = run(gen, "--src", src, "--dst", inside, "--ids", "c", cwd=tmp_path)
    assert r.returncode == 2 and "refusing" in r.stdout
    assert not inside.exists()
    # a sibling destination gets five frequencies with the documented signatures
    dst = tmp_path / "tree" / "mock"
    r = run(gen, "--src", src, "--dst", dst, "--ids", "c", cwd=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "wrote 20 CSV files and 5 sidecars to " in r.stdout
    base = str(dst / "waveguides")
    grid = hfss.discover_frequencies("c", base)
    assert [f for f, _ in grid] == [50.0, 150.0, 500.0, 1500.0, 5000.0]
    for (f, stem), scale in zip(grid, (0.2, 0.4, 1.0, 0.6, 0.8)):
        ds = hfss.load_dataset(stem, base)[(0.0, 180.0)]
        assert ds.T0 == pytest.approx(0.5 * scale) and ds.T1 == pytest.approx(0.25 * scale)
        back = ds.theta_deg > 90.
        amps = np.abs(np.concatenate([ds.F0_theta[back], ds.F0_phi[back], ds.F1_theta[back], ds.F1_phi[back]]))
        truncated = f in (150.0, 1500.0)
        assert np.all(amps == 0) if truncated else np.all(amps > 0)
    for f, stem in grid:
        sc = sidecar.check_full(base, "c", stem, f)      # every mock frequency has a valid sidecar
        assert sc["frequency_ghz"] == f and sc["provenance"]["mock_of"] == "c_500GHz"
        assert sc["modes"]["propagating_count"] == sidecar.modes_for(sc["exit_field"]["cross_section"], f)["propagating_count"]


def test_mock_generator_needs_source_sidecar(repo_root, tmp_path):
    # A source dataset without its sidecar: rc 2, the BBR024 message, nothing written.
    gen = repo_root / "validation/Scripts/make_mock_hfss_frequencies.py"
    src = tmp_path / "tree" / "data" / "waveguides"
    write_source(src, "c")
    missing = src / "c_500GHz.dataset.json"
    missing.unlink()
    dst = tmp_path / "tree" / "mock"
    r = run(gen, "--src", src, "--dst", dst, "--ids", "c", cwd=tmp_path)
    assert r.returncode == 2, r.stdout + r.stderr
    assert f"source sidecar: BBR024: {missing.resolve()}: missing" in r.stdout
    assert not dst.exists()


def test_mock_round_gap(repo_root, tmp_path, data_root):
    gen = repo_root / "validation/Scripts/make_mock_round_gap.py"
    real = os.path.join(data_root, "waveguides")
    r = run(gen, "--real", real, "--dst", tmp_path / "rg", cwd=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    base = str(tmp_path / "rg" / "waveguides")
    sc = sidecar.check_full(base, "RoundGap_r50um", "RoundGap_r50um_2000GHz", 2000.0)
    assert sc["modes"]["mode"] == "TE11" and sc["modes"]["propagating_count"] == 1
    assert sc["exit_field"]["outside_points"] == "omitted" and sc["symmetry"]["rotational"] is True
    ds = hfss.load_dataset("RoundGap_r50um_2000GHz", base)
    assert ds[(0.0, 180.0)].T0 == pytest.approx(0.8) and ds[(0.0, 180.0)].T1 == pytest.approx(0.6)
    assert abs(ds[(45.0, 180.0)].rho.real) > 0.1 and ds[(0.0, 180.0)].rho == 0
    import pandas as pd
    wg = pd.read_csv(os.path.join(base, "RoundGap_r50um_2000GHz_Ephi=1", "waveguide.csv"))
    assert np.hypot(wg.Y, wg.Z).max() == pytest.approx(5e-5, rel=1e-12)    # rim points kept
    rim = np.isclose(np.hypot(wg.Y, wg.Z), 5e-5, rtol=1e-12, atol=0.0)
    per_key = wg.assign(rim=rim).groupby(["IWavePhi", "IWaveTheta"]).rim.sum()
    assert (per_key == 20).all() and sc["exit_field"]["points_per_key_retained"] == 1961
    assert os.path.islink(os.path.join(base, "InfParallelPlate_crack1Rohan_500GHz_Ephi=0"))
    assert os.path.isfile(sidecar.path_for(base, "InfParallelPlate_crack2_500GHz"))
    r = run(gen, "--real", real, "--dst", os.path.join(data_root, "rg"), cwd=tmp_path)
    assert r.returncode == 2 and "refusing" in r.stdout
