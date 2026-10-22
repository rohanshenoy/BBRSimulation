"""bbrsim.hfss on synthetic few-row dataset trees (and one real-data key)."""
import os

import numpy as np
import pandas as pd
import pytest

from bbrsim import hfss

FF_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag\n"
WG_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag\n"
KEYS = [(0.0, 180.0), (0.0, 135.0)]
GRID = [(-90.0, 90.0), (0.0, 90.0), (45.0, 60.0), (0.0, 120.0)]   # (Phi, Theta)


def write_dataset(base, stem, T=(0.5, 0.9), ratio_over_one=False, freq="500GHz", grid1=None,
                  fields=((0, 1, 0), (0, 0, 1))):
    """Two keys, a 4-point far field and one exit point per key; fields[e] is the
    complex (Ex, Ey, Ez) of the Ephi=e exit point (default orthogonal: rho = 0; a
    transmitting polarization needs a nonzero one, BBR027)."""
    for e in (0, 1):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d)
        g = GRID if (e == 0 or grid1 is None) else grid1
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(FF_HDR)
            for ip, it in KEYS:
                for p, t in g:
                    fh.write(f"{freq},{e},{ip},{it},{p},{t},{0.1*(e+1)},0.0,{0.2},{0.05}\n")
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(WG_HDR)
            f = [complex(c) for c in fields[e]]
            for ip, it in KEYS:
                out = 1.0545 if ratio_over_one else T[e]
                fh.write(f"{freq},{e},{ip},{it},{out},1.0,0,0,0,"
                         + ",".join(str(c.real) for c in f) + ","
                         + ",".join(str(c.imag) for c in f) + "\n")


def edit_last_row(path, column, value):
    """Set one column of the last data row of a CSV (the fixtures' second key)."""
    lines = path.read_text().splitlines()
    v = lines[-1].split(",")
    v[column] = value
    lines[-1] = ",".join(v)
    path.write_text("\n".join(lines) + "\n")


def lam_max(T0, T1, rho):
    """Largest coherent T over linear polarizations: the top eigenvalue of
    [[T0, c], [c, T1]] with c = sqrt(T0 T1) Re rho."""
    c = np.sqrt(T0 * T1) * rho.real
    return (T0 + T1) / 2 + np.sqrt(((T0 - T1) / 2) ** 2 + c ** 2)


def raw_ratios(base, stem):
    """{key: (T0, T1)}, OutgoingPower/IngoingPower of each key's first waveguide row per
    Ephi, read from the CSVs (not through load_dataset)."""
    out = {}
    for e in (0, 1):
        wg = pd.read_csv(os.path.join(base, f"{stem}_Ephi={e}", "waveguide.csv"),
                         usecols=["IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower"])
        for (p, t), g in wg.groupby(["IWavePhi", "IWaveTheta"], sort=False):
            r = g.iloc[0]
            key = (round(float(p), hfss.KEY_ROUND), round(float(t), hfss.KEY_ROUND))
            out.setdefault(key, [0.0, 0.0])[e] = (
                float(r.OutgoingPower) / float(r.IngoingPower) if r.IngoingPower > 0 else 0.0)
    return {k: tuple(v) for k, v in out.items()}


# --- frequency discovery / selection ------------------------------------------
def test_discover_sorted_verbatim_token(tmp_path):
    for tok in ("500", "1.5e3", "50"):
        write_dataset(tmp_path, f"crack_{tok}GHz")
    (tmp_path / "crack_notanumberGHz_Ephi=0").mkdir()
    (tmp_path / "crack_0GHz_Ephi=0").mkdir()
    (tmp_path / "crack_700GHz_Ephi=0").write_text("a file, not a dir")
    assert hfss.discover_frequencies("crack", str(tmp_path)) == [
        (50.0, "crack_50GHz"), (500.0, "crack_500GHz"), (1500.0, "crack_1.5e3GHz")]


def test_discover_legacy_dir_raises(tmp_path):
    (tmp_path / "crack_Ephi=0").mkdir()
    with pytest.raises(ValueError, match="legacy"):
        hfss.discover_frequencies("crack", str(tmp_path))


def test_discover_duplicate_raises(tmp_path):
    (tmp_path / "crack_1500GHz_Ephi=0").mkdir()
    (tmp_path / "crack_1.5e3GHz_Ephi=0").mkdir()
    with pytest.raises(ValueError, match="duplicate"):
        hfss.discover_frequencies("crack", str(tmp_path))


def test_discover_nothing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        hfss.discover_frequencies("crack", str(tmp_path))


def test_discover_prefix_is_not_a_substring_match(tmp_path):
    (tmp_path / "crack2_500GHz_Ephi=0").mkdir()
    with pytest.raises(FileNotFoundError):
        hfss.discover_frequencies("crack", str(tmp_path))


GRID5 = [(50.0, "a"), (150.0, "b"), (500.0, "c"), (1500.0, "d"), (5000.0, "e")]


@pytest.mark.parametrize("nu,expected", [
    (20.0, (50.0, "a", -1)), (10000.0, (5000.0, "e", +1)),
    (50.0, (50.0, "a", 0)), (5000.0, (5000.0, "e", 0)),
    (np.sqrt(50 * 150), (50.0, "a", 0)),            # exact log tie -> lower
    (np.sqrt(50 * 150) * 1.01, (150.0, "b", 0)),
])
def test_select_frequency(nu, expected):
    assert hfss.select_frequency(GRID5, nu) == expected


def test_single_grid_never_clamps():
    assert hfss.select_frequency([(500.0, "c")], 1e6) == (500.0, "c", 0)


@pytest.mark.parametrize("lo,hi", [(50.0, 150.0), (150.0, 500.0), (500.0, 1500.0), (1500.0, 5000.0)])
def test_select_frequency_log_midpoints(lo, hi):
    # sqrt(lo*hi) is an exact tie in log10 distance in binary for all four midpoints:
    # it goes to the lower frequency; 1 % either side picks the nearer one.
    mid = np.sqrt(lo * hi)
    assert hfss.select_frequency(GRID5, mid)[0] == lo
    assert hfss.select_frequency(GRID5, mid * 0.99)[0] == lo
    assert hfss.select_frequency(GRID5, mid * 1.01)[0] == hi


def test_photon_frequency_500GHz():
    assert hfss.photon_frequency_GHz(2.067834e-3) == pytest.approx(500.0, abs=1e-3)


# Directory-name tokens as C strtod full-parses them (BBRCrackLibrary::Discover):
# verified against the macOS libc strtod.
@pytest.mark.parametrize("token,value", [
    ("500", 500.0), ("1.5e3", 1500.0), (" 500", 500.0), ("\t500", 500.0), ("+500", 500.0),
    (".5e3", 500.0), ("1.", 1.0), ("0x1f4", 500.0), ("+0x1f4", 500.0), ("0x1.f4p8", 500.0),
    ("inf", np.inf), ("INF", np.inf), ("Infinity", np.inf),
    ("500 ", None), ("5_00", None), ("", None), ("-5", None), ("0", None), ("nan", None),
    ("nan(1)", None), ("1e", None), ("0x", None), ("0x1p", None), ("-inf", None),
    ("1e-999", None), ("500\n", None), ("٥٠٠", None),
])
def test_token_corpus_matches_cpp(token, value):
    got = hfss._parse_stem_token(f"c_{token}GHz_Ephi=0", "c_", "GHz_Ephi=0")
    assert (got is None) if value is None else (got == (value, token))


def test_photon_frequency_scalar_and_array():
    assert isinstance(hfss.photon_frequency_GHz(2.067834e-3), float)
    assert hfss.photon_frequency_GHz(np.array([2.067834e-3])).shape == (1,)


def test_photon_frequency_matches_clhep_arithmetic():
    # CLHEP builds h_Planck = 6.62607015e-34 * joule * s in MeV ns, with joule =
    # 1e-6 MeV / e_SI and s = 1e9 ns; BBSimOpBoundaryProcess computes
    # nu = E[MeV] / h_Planck / (1e9 * hertz), and 1e9 * hertz is exactly 1. The
    # literals below restate that arithmetic independently of the module.
    assert hfss.H_PLANCK_MEV_NS == 4.135667696923858e-12
    assert hfss.H_EV_S == 4.135667696923858e-15
    rng = np.random.default_rng(7)
    e_MeV = rng.uniform(4.14e-5, 8.27e-2, 10000) * 1e-6     # the emitter's band, as Geant4 holds it
    recorded_eV = e_MeV / 1e-6                               # the ntuple stores E / eV
    nu_cpp = e_MeV / 4.135667696923858e-12 / 1.0
    assert np.array_equal(hfss.photon_frequency_GHz(recorded_eV), nu_cpp)


def test_discover_datasets_lists_every_id(tmp_path):
    for name in ("InfParallelPlate_crack1Rohan_180GHz", "InfParallelPlate_crack1Rohan_1.01e3GHz",
                 "InfParallelPlate_crack2_500GHz"):
        (tmp_path / f"{name}_Ephi=0").mkdir()
        (tmp_path / f"{name}_Ephi=1").mkdir()
        (tmp_path / f"{name}.dataset.json").write_text("{}")
    (tmp_path / "c_notanumberGHz_Ephi=0").mkdir()            # skipped, as BBRCrackLibrary skips it
    (tmp_path / "c_700GHz_Ephi=0").write_text("a file, not a dir")
    assert hfss.discover_datasets(str(tmp_path)) == {
        "InfParallelPlate_crack1Rohan": [(180.0, "InfParallelPlate_crack1Rohan_180GHz"),
                                         (1010.0, "InfParallelPlate_crack1Rohan_1.01e3GHz")],
        "InfParallelPlate_crack2": [(500.0, "InfParallelPlate_crack2_500GHz")]}


def test_discover_datasets_legacy_dir_raises(tmp_path):
    (tmp_path / "InfParallelPlate_crack2_Ephi=0").mkdir()
    with pytest.raises(ValueError, match="legacy"):
        hfss.discover_datasets(str(tmp_path))


def test_discover_datasets_duplicate_raises(tmp_path):
    (tmp_path / "c_1500GHz_Ephi=0").mkdir()
    (tmp_path / "c_1.5e3GHz_Ephi=0").mkdir()
    with pytest.raises(ValueError, match="duplicate"):
        hfss.discover_datasets(str(tmp_path))


def test_discover_datasets_empty_tree(tmp_path):
    assert hfss.discover_datasets(str(tmp_path)) == {}


# --- dataset load ----------------------------------------------------------------
def test_load_normalizes_ratio_above_one(tmp_path):
    # A key whose lam_max exceeds 1 has T0 and T1 divided by it (the BBRHFSSData
    # constructor). rho = 0 and T0 = T1 = 1.0545: lam_max = 1.0545, both become 1.
    write_dataset(tmp_path / "cap", "c_500GHz", ratio_over_one=True)
    ds = hfss.load_dataset("c_500GHz", str(tmp_path / "cap"))
    assert set(ds) == {(0.0, 180.0), (0.0, 135.0)}
    assert all(d.T0 == 1.0 and d.T1 == 1.0 for d in ds.values())
    # rho = 0 with T0 = 1.25, T1 = 0.5: T0 -> 1 and T1 -> 0.4 (the per-polarization cap
    # this replaces left T1 at 0.5).
    write_dataset(tmp_path / "one", "c_500GHz", T=(1.25, 0.5))
    d = hfss.load_dataset("c_500GHz", str(tmp_path / "one"))[(0.0, 180.0)]
    assert (d.T0, d.T1) == (1.0, 0.4)
    # rho = 1 with T0 = 0.8, T1 = 0.4: the in-phase T would be 1.2 although neither T0 nor
    # T1 exceeds 1; after the division it is 1 and T1/T0 is kept.
    write_dataset(tmp_path / "sum", "c_500GHz", T=(0.8, 0.4), fields=((0, 1, 0), (0, 1, 0)))
    d = hfss.load_dataset("c_500GHz", str(tmp_path / "sum"))[(0.0, 180.0)]
    assert d.T0 == pytest.approx(0.8 / 1.2, abs=1e-15) and d.T1 / d.T0 == pytest.approx(0.5, abs=1e-15)
    assert hfss.transmittance(d, (2 / 3) ** .5, (1 / 3) ** .5) == pytest.approx(1.0, abs=1e-15)
    assert hfss.transmittance(d, (1 / 3) ** .5, -(2 / 3) ** .5) == pytest.approx(0.0, abs=1e-15)
    # lam_max <= 1 (rho = 1, 0.5 + 0.4 = 0.9): untouched, bit for bit.
    write_dataset(tmp_path / "low", "c_500GHz", T=(0.5, 0.4), fields=((0, 1, 0), (0, 1, 0)))
    d = hfss.load_dataset("c_500GHz", str(tmp_path / "low"))[(0.0, 180.0)]
    assert (d.T0, d.T1) == (0.5, 0.4)


def test_load_rejects_mismatched_grids(tmp_path):
    write_dataset(tmp_path, "c_500GHz", grid1=[(-90.0, 90.0), (0.0, 90.0), (45.0, 61.0), (0.0, 120.0)])
    with pytest.raises(ValueError, match="grids differ"):
        hfss.load_dataset("c_500GHz", str(tmp_path))


def test_nearest_key():
    ds = {(0.0, 180.0): None, (0.0, 135.0): None}
    assert hfss.nearest_key(ds, 10.0, 140.0) == (0.0, 135.0)


def test_nearest_key_tie_is_map_order():
    # (45, 180) is equidistant from both keys; std::map order (ascending) scans
    # (0, 180) first and keeps it, whatever order the CSV listed the keys in.
    ds = {(90.0, 180.0): None, (0.0, 180.0): None}
    assert hfss.nearest_key(ds, 45.0, 180.0) == (0.0, 180.0)


@pytest.mark.parametrize("name", ["far_field.csv", "waveguide.csv"])
def test_load_checks_every_freq_row(tmp_path, name):
    write_dataset(tmp_path, "c_500GHz")
    edit_last_row(tmp_path / "c_500GHz_Ephi=1" / name, 0, "140GHz")
    with pytest.raises(ValueError, match="BBR009"):
        hfss.load_dataset("c_500GHz", str(tmp_path))


def test_load_accepts_freq_units(tmp_path):
    for i, tok in enumerate(("0.5THz", "500000MHz", " 500 ghz")):
        write_dataset(tmp_path / str(i), "c_500GHz", freq=tok)
        assert set(hfss.load_dataset("c_500GHz", str(tmp_path / str(i)))) == set(KEYS)


def test_load_rejects_exit_point_mismatch(tmp_path):
    write_dataset(tmp_path / "y", "c_500GHz")
    edit_last_row(tmp_path / "y" / "c_500GHz_Ephi=1" / "waveguide.csv", 7, "0.001")   # Y
    with pytest.raises(ValueError, match="BBR012"):
        hfss.load_dataset("c_500GHz", str(tmp_path / "y"))
    write_dataset(tmp_path / "n", "c_500GHz")
    wg = tmp_path / "n" / "c_500GHz_Ephi=1" / "waveguide.csv"
    wg.write_text(wg.read_text() + wg.read_text().splitlines()[-1] + "\n")   # one Ephi=1 row too many
    with pytest.raises(ValueError, match="BBR012"):
        hfss.load_dataset("c_500GHz", str(tmp_path / "n"))


def test_load_rejects_key_mismatch(tmp_path):
    # waveguide.csv at a key far_field.csv lacks: BBR007, as in the C++ (not an IndexError).
    write_dataset(tmp_path / "k", "c_500GHz")
    for e in (0, 1):
        edit_last_row(tmp_path / "k" / f"c_500GHz_Ephi={e}" / "waveguide.csv", 3, "90")
    with pytest.raises(ValueError, match="BBR007"):
        hfss.load_dataset("c_500GHz", str(tmp_path / "k"))
    # an Ephi=1 far-field row at a key with no Ephi=0 rows: BBR012
    write_dataset(tmp_path / "x", "c_500GHz")
    ff1 = tmp_path / "x" / "c_500GHz_Ephi=1" / "far_field.csv"
    ff1.write_text(ff1.read_text() + "500GHz,1,0,90,0,90,0.2,0.0,0.2,0.05\n")
    with pytest.raises(ValueError, match="BBR012"):
        hfss.load_dataset("c_500GHz", str(tmp_path / "x"))


def test_load_rejects_nan_field(tmp_path):
    # A field that is not a finite number is BBR013, as the C++ Num raises it.
    for i, v in enumerate(("nan", "abc", "inf")):
        write_dataset(tmp_path / str(i), "c_500GHz")
        edit_last_row(tmp_path / str(i) / "c_500GHz_Ephi=0" / "far_field.csv", 8, v)   # rEtheta_real
        with pytest.raises(ValueError, match="BBR013"):
            hfss.load_dataset("c_500GHz", str(tmp_path / str(i)))


# --- geometry: fold, basis, weights ------------------------------------------------
def test_fold_normal_incidence():
    inc = hfss.fold_incidence((1, 0, 0))
    assert inc.theta_deg == pytest.approx(180.0) and inc.phi_deg == 0.0 and (inc.sy, inc.sx) == (1., 1.)


def test_fold_signed_zero_snap():
    inc = hfss.fold_incidence((1, -0.0, -0.5))
    assert inc.phi_raw_deg == pytest.approx(180.0) and inc.sy == 1.0 and inc.sx == -1.0


def test_fold_backward_photon_flips_normal():
    inc = hfss.fold_incidence((-1, 0, 0))
    assert np.allclose(inc.normal, (-1, 0, 0)) and inc.theta_deg == pytest.approx(180.0)


def test_incoming_basis_is_orthonormal_and_transverse():
    # Pins the deliberate +sin(theta) sign (a documented convention): with the other sign,
    # e_theta . k != 0 for any oblique photon.
    rng = np.random.default_rng(1)
    for _ in range(500):
        k = rng.normal(size=3)
        k[0] = abs(k[0]) + 1e-3
        k /= np.linalg.norm(k)
        inc = hfss.fold_incidence(k)
        assert 0.0 <= inc.phi_deg <= 90.0
        et, ep = hfss.incoming_basis(inc)
        assert abs(np.dot(et, k)) < 1e-12 and abs(np.dot(ep, k)) < 1e-12
        assert abs(np.dot(et, ep)) < 1e-12
        assert np.linalg.norm(et) == pytest.approx(1) and np.linalg.norm(ep) == pytest.approx(1)


def test_polarization_components_fallback():
    inc = hfss.fold_incidence((1, 0, 0))
    assert hfss.polarization_components(np.array([1.0, 0, 0]), inc) == pytest.approx((2 ** -0.5, 2 ** -0.5))


def test_transmittance_clamp():
    class D: T0, T1, rho = 1.0, 1.0, 0j
    assert hfss.transmittance(D, 1.0, 1.0) == 1.0


R2 = 2 ** -0.5


def test_transmittance_cross_term(tmp_path):
    # rho = 1 (identical exit fields): the transmitted wave is the projection onto one
    # axis, T = (Et sqrt(T0) + Ep sqrt(T1))^2, so with T0 = T1 the in-phase diagonal
    # gives T0 + T1 and the orthogonal one 0.
    write_dataset(tmp_path / "eq", "c_500GHz", T=(0.4, 0.4), fields=((0, 1, 0), (0, 1, 0)))
    d = hfss.load_dataset("c_500GHz", str(tmp_path / "eq"))[(0.0, 180.0)]
    assert d.rho == pytest.approx(1.0)
    assert hfss.transmittance(d, R2, R2) == pytest.approx(0.8, abs=1e-12)
    assert hfss.transmittance(d, R2, -R2) == pytest.approx(0.0, abs=1e-12)
    write_dataset(tmp_path / "ne", "c_500GHz", T=(0.3, 0.5), fields=((0, 1, 0), (0, 1, 0)))
    d = hfss.load_dataset("c_500GHz", str(tmp_path / "ne"))[(0.0, 180.0)]
    assert hfss.transmittance(d, R2, R2) == pytest.approx((0.3 ** .5 + 0.5 ** .5) ** 2 / 2, abs=1e-12)
    assert hfss.transmittance(d, R2, -R2) == pytest.approx((0.3 ** .5 - 0.5 ** .5) ** 2 / 2, abs=1e-12)
    # rho = 0 (orthogonal fields) and rho = i (quadrature): no cross term.
    for i, f1 in enumerate(((0, 0, 1), (0, 1j, 0))):
        write_dataset(tmp_path / f"o{i}", "c_500GHz", T=(0.4, 0.4), fields=((0, 1, 0), f1))
        d = hfss.load_dataset("c_500GHz", str(tmp_path / f"o{i}"))[(0.0, 180.0)]
        assert d.rho.real == pytest.approx(0.0, abs=1e-15)
        assert hfss.transmittance(d, R2, R2) == pytest.approx(0.4, abs=1e-12)
        assert hfss.transmittance(d, R2, -R2) == pytest.approx(0.4, abs=1e-12)


@pytest.mark.parametrize("f1", [(0, 1, 0), (0, 0, 1), (0, -1, 0), (0, 0.6, 0.8j)])
def test_random_polarization_mixture_unchanged_by_cross_term(tmp_path, f1):
    # The psi-average of 2 Et Ep sqrt(T0 T1) Re rho is zero, so the unpolarized T is
    # (T0 + T1)/2 for any rho, as long as the clamp at 1 never acts (T0 + T1 <= 1).
    write_dataset(tmp_path, "c_500GHz", T=(0.3, 0.5), fields=((0, 1, 0), f1))
    d = hfss.load_dataset("c_500GHz", str(tmp_path))[(0.0, 180.0)]
    assert hfss.random_polarization_mixture(d)[0] == pytest.approx(0.4, abs=1e-12)


def test_real_data_diagonal_key(data_root):
    # (45, 135) of crack1: the two basis exit fields are 99.96 % correlated, so an
    # in-phase diagonal polarization transmits ~T0 + T1 and the orthogonal one ~0.
    base = os.path.join(data_root, "waveguides")
    # The random-polarization mean at (45, 180), where the raw in-phase T exceeds 1 (T0 + T1
    # = 1.025 / 1.006, Re rho ~ 1) while neither T0 nor T1 does: the load-time normalization
    # moves it from (T0 + T1)/2 = 0.512645 / 0.502883 to 0.500094 / 0.500033.
    mean_45_180 = {"InfParallelPlate_crack1Rohan_500GHz": 0.500094, "InfParallelPlate_crack2_500GHz": 0.500033}
    # The keys whose raw lam_max exceeds 1: (45, 180), and those where T0 or T1 alone does.
    normalized = {"InfParallelPlate_crack1Rohan_500GHz": {(0.0, 180.0), (45.0, 180.0)},
                  "InfParallelPlate_crack2_500GHz": {(0.0, 180.0), (45.0, 180.0), (90.0, 180.0)}}
    for stem in mean_45_180:
        if not os.path.isdir(os.path.join(base, f"{stem}_Ephi=0")):
            pytest.skip(f"no {stem} dataset under {base}")
    d = hfss.load_dataset("InfParallelPlate_crack1Rohan_500GHz", base)[(45.0, 135.0)]
    assert d.T0 == pytest.approx(0.4073, abs=5e-4) and d.T1 == pytest.approx(0.3939, abs=5e-4)
    assert hfss.transmittance(d, R2, R2) == pytest.approx(d.T0 + d.T1, abs=0.002)
    assert hfss.transmittance(d, R2, -R2) == pytest.approx(0.0, abs=0.001)
    # Every key of both cracks (spec P1). At IWavePhi 0 and 90 the cross term stays within
    # 2 |Et Ep| sqrt(T0 T1) (at most 4e-5 here). After the load no key's lam_max exceeds 1,
    # so the clamp never acts and the random-polarization mean is (T0 + T1)/2 at every key.
    # Against the raw first-row ratios (read here from the CSVs): a key with raw lam_max <= 1
    # is untouched bit for bit; one where T0 or T1 alone exceeds 1 matches the old
    # per-polarization cap to 1e-9; only (45, 180) moves beyond that.
    for stem, mean45 in mean_45_180.items():
        ds = hfss.load_dataset(stem, base)
        raw = raw_ratios(base, stem)
        assert len(ds) == 15 and set(raw) == set(ds)
        hit = set()
        for key, d in ds.items():
            if key[0] in (0.0, 90.0):
                for Et, Ep in ((R2, R2), (R2, -R2)):
                    cross = hfss.transmittance(d, Et, Ep) - (Et ** 2 * d.T0 + Ep ** 2 * d.T1)
                    assert abs(cross) <= 2 * abs(Et * Ep) * np.sqrt(d.T0 * d.T1) + 1e-15, (stem, key)
            assert lam_max(d.T0, d.T1, d.rho) <= 1 + 1e-12, (stem, key)
            mean = hfss.random_polarization_mixture(d)[0]
            assert mean == pytest.approx((d.T0 + d.T1) / 2, abs=1e-12), (stem, key)
            r0, r1 = raw[key]
            if lam_max(r0, r1, d.rho) <= 1:
                assert (d.T0, d.T1) == (r0, r1), (stem, key)
                continue
            hit.add(key)
            if key == (45.0, 180.0):
                assert mean == pytest.approx(mean45, abs=1e-6), (stem, key)
            else:
                assert abs(d.T0 - min(1., r0)) <= 1e-9 and abs(d.T1 - min(1., r1)) <= 1e-9, (stem, key)
        assert hit == normalized[stem]
        assert ds[(0.0, 180.0)].T0 == pytest.approx(1.0, abs=1e-9)


def test_direction_weights(tmp_path):
    write_dataset(tmp_path, "c_500GHz")
    ds = hfss.load_dataset("c_500GHz", str(tmp_path))[(0.0, 180.0)]
    w = hfss.direction_weights(ds, 1.0, 0.0)
    assert w.sum() == pytest.approx(1.0)
    assert w[0] == 0.0                     # Phi = -90: in the exit-face plane
    inc = hfss.fold_incidence((1, 0, 0))
    dirs = hfss.outgoing_directions(ds, inc)
    assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0)
    assert np.all(dirs[w > 0] @ inc.normal > 0)


def write_exit_dataset(base, stem, points, f0, f1):
    """Both KEYS with the 4-point far field of write_dataset and the given exit points
    (Y, Z in metres, X = 0); f0[j] / f1[j] are the complex (Ex, Ey, Ez) at point j."""
    for e, fields in ((0, f0), (1, f1)):
        d = os.path.join(base, f"{stem}_Ephi={e}")
        os.makedirs(d)
        with open(os.path.join(d, "far_field.csv"), "w") as fh:
            fh.write(FF_HDR)
            for ip, it in KEYS:
                for p, t in GRID:
                    fh.write(f"500GHz,{e},{ip},{it},{p},{t},0.1,0.0,0.2,0.05\n")
        with open(os.path.join(d, "waveguide.csv"), "w") as fh:
            fh.write(WG_HDR)
            for ip, it in KEYS:
                for (y, z), f in zip(points, fields):
                    f = [complex(c) for c in f]
                    fh.write(f"500GHz,{e},{ip},{it},0.4,1.0,0,{y},{z},"
                             + ",".join(str(c.real) for c in f) + ","
                             + ",".join(str(c.imag) for c in f) + "\n")


def test_exit_position_weights(tmp_path):
    # Three exit points; weight_j = |Et E0_j + Ep E1_j|^2 summed over components
    # (BBRHFSSData::SampleExitPosition), normalized; BBR027 when every weight is zero.
    pts = [(-1e-5, 0.0), (0.0, 2e-5), (3e-5, -4e-5)]
    f0 = [(0, 1, 0), (0, 0, 1), (0, 0, 0)]
    f1 = [(0, 1, 0), (0, 0, -1), (0, 0, 2j)]
    write_exit_dataset(tmp_path, "c_500GHz", pts, f0, f1)
    d = hfss.load_dataset("c_500GHz", str(tmp_path))[(0.0, 135.0)]
    assert np.allclose(d.exit_y_m, [p[0] for p in pts]) and np.allclose(d.exit_z_m, [p[1] for p in pts])
    assert d.E0.shape == (3, 3) and d.E1[2, 2] == 2j
    assert np.allclose(hfss.exit_position_weights(d, 1.0, 0.0), [0.5, 0.5, 0.0])
    assert np.allclose(hfss.exit_position_weights(d, 0.0, 1.0), [1 / 6, 1 / 6, 4 / 6])
    # The diagonal adds the fields coherently: point 0 in phase (|2 R2|^2 = 2), point 1
    # cancels (0), point 2 in quadrature (|2j R2|^2 = 2).
    assert np.allclose(hfss.exit_position_weights(d, R2, R2), [0.5, 0.0, 0.5])
    assert np.allclose(hfss.exit_position_weights(d, R2, -R2), [0.0, 0.5, 0.5])
    # CDX-008: T = 0.4 > 0 with no exit field at all is refused at load, not sampled uniformly.
    write_exit_dataset(tmp_path / "zero", "c_500GHz", pts, [(0, 0, 0)] * 3, [(0, 0, 0)] * 3)
    with pytest.raises(ValueError, match=r"^BBR027: c_500GHz key \(0.0, 180.0\): T0 = 0.4 > 0, but its exit field "
                                         r"is zero at every exit point"):
        hfss.load_dataset("c_500GHz", str(tmp_path / "zero"))


def test_weights_refuse_cancelling_fields(tmp_path):
    """CDX-008: the basis fields cancel for one mixed polarization; the mirror raises BBR027, as the
    samplers do, instead of returning uniform weights."""
    write_exit_dataset(tmp_path, "c_500GHz", [(0.0, 0.0)], [(0, 1, 0)], [(0, -1, 0)])
    d = hfss.load_dataset("c_500GHz", str(tmp_path))[(0.0, 180.0)]
    with pytest.raises(ValueError, match=r"^BBR027: no exit field for \(E_theta, E_phi\)"):
        hfss.exit_position_weights(d, R2, R2)
    with pytest.raises(ValueError, match=r"^BBR027: no outward far-field intensity"):
        hfss.direction_weights(d, R2, -R2)   # the two polarizations' far fields are equal: they cancel here
    assert np.allclose(hfss.exit_position_weights(d, 1.0, 0.0), [1.0])


def test_random_polarization_mean_is_average(tmp_path):
    write_dataset(tmp_path, "c_500GHz", T=(0.3, 0.8))
    ds = hfss.load_dataset("c_500GHz", str(tmp_path))[(0.0, 180.0)]
    T_mean, W = hfss.random_polarization_mixture(ds)
    assert T_mean == pytest.approx(0.55, abs=1e-12) and W.sum() == pytest.approx(1.0)


def test_reflected_direction():
    r = hfss.reflected_direction((0.6, 0.8, 0.0))
    assert np.allclose(r, (-0.6, 0.8, 0.0))


def test_binned_expectation_conserves_weight():
    out = hfss.binned_expectation(np.array([-5.0, 0.5, 1.5, 99.0]), np.ones(4), np.array([0.0, 1.0, 2.0]))
    assert out.tolist() == [2.0, 2.0]


# --- schema 2.x (Task 6): the mirror of BBRHFSSData's 2.x branch ---------------------------------------
# Fixture: schema2_fixture.write_nine_key_grid2, the Python twin of the C++ NineKeyGrid2.
import schema2_fixture as fx   # noqa: E402

from bbrsim import sidecar   # noqa: E402

NODE_W = (0.5, 1.0, 0.5)   # the trapezoid factor of a 3-node axis, by grid index


def s2_load(tmp_path, tflux=None, edit=None):
    """Write the nine-key dataset and its 2.0 sidecar under tmp_path; return its AngleDataset map."""
    base = str(tmp_path)
    stem = fx.write_nine_key_grid2(base)
    if edit:
        edit(base, stem)
    sidecar.write(base, stem, fx.sidecar2(tflux=tflux))
    return hfss.load_dataset(stem, base)


def s2_fields(e, i):
    """(Ex, Ey, Ez) complex of exit point i (row-major over y then z), polarization e, as the fixture writes
    them: six values (i+1+f) * 1e-3 (Ephi=0) or (i+6-f) * 1e-3 (Ephi=1) in the 2.0 order Re/Im pairs."""
    v = [(i + 1 + f if e == 0 else i + 6 - f) * 1e-3 for f in range(6)]
    return np.array([v[0] + 1j * v[1], v[2] + 1j * v[3], v[4] + 1j * v[5]])


def s2_w(i):
    return NODE_W[i // 3] * NODE_W[i % 3]


def test_s2_positions(tmp_path):
    ds = s2_load(tmp_path)
    d = ds[(45.0, 135.0)]
    assert d.schema == 2 and not d.dark0 and not d.dark1
    assert np.array_equal(d.w, [s2_w(i) for i in range(9)])
    # E_theta only: weight_i = w_i (v9^2 + v10^2 + v11^2 + v12^2), the four values after X, Y, Z that the 2.0
    # map names Ex_real, Ex_imag, Ey_real, Ey_imag (the field-list positions 9 to 12).
    want = np.array([s2_w(i) * sum(((i + 1 + f) * 1e-3) ** 2 for f in range(4)) for i in range(9)])
    got = hfss.exit_position_weights(d, 1.0, 0.0)
    assert got == pytest.approx(want / want.sum(), rel=1e-12, abs=0)
    # the 1.0 map on the same columns (Ex_re, Ey_re, Ez_re, Ex_im, Ey_im, Ez_im at 9..14, uniform weights, all
    # three components) gives a different distribution
    old = np.array([sum(((i + 1 + f) * 1e-3) ** 2 for f in range(6)) for i in range(9)])
    assert not np.allclose(got, old / old.sum(), rtol=1e-3, atol=0)
    # a mixed polarization: w_i |cT Et0 + cP Et1|^2 over (Ex, Ey)
    mix = np.array([s2_w(i) * np.sum(np.abs(R2 * s2_fields(0, i)[:2] + R2 * s2_fields(1, i)[:2]) ** 2) for i in range(9)])
    assert hfss.exit_position_weights(d, R2, R2) == pytest.approx(mix / mix.sum(), rel=1e-12, abs=0)
    # Ez is not read: the exit sampler never sees the third component
    assert np.allclose(d.E0[:, 2], [s2_fields(0, i)[2] for i in range(9)])


def test_s2_transmittance_and_cap(tmp_path):
    ds = s2_load(tmp_path)
    ck90 = sidecar.cell_90([90.0, 135.0, 180.0])
    d = ds[(0.0, 90.0)]
    # T_p = OutgoingPower / (IngoingPower * InterceptFactor) from the key's first row, kept above 1
    assert d.T0 == pytest.approx(0.6 / ck90, rel=1e-12) and d.T0 > 1
    assert d.T1 == pytest.approx(0.3 / ck90, rel=1e-12) and d.T1 > 1
    assert hfss.transmittance(d, 1.0, 0.0) == 1.0          # the per-photon cap is the clamp
    assert hfss.transmittance(d, 0.0, 1.0) == 1.0
    d135 = ds[(0.0, 135.0)]
    assert d135.T0 == pytest.approx(0.6 / np.cos(np.pi / 4), rel=1e-12)
    assert hfss.transmittance(d135, 1.0, 0.0) == pytest.approx(min(1.0, d135.T0), rel=1e-12)
    d180 = ds[(45.0, 180.0)]
    assert (d180.T0, d180.T1) == pytest.approx((0.6, 0.3), rel=1e-12)


def test_s2_no_rescale(tmp_path):
    ds = s2_load(tmp_path)
    for (phi, theta), d in ds.items():
        ck = sidecar.intercept_factor(theta, [90.0, 135.0, 180.0])
        assert d.T0 == pytest.approx(0.6 / ck, rel=1e-12), (phi, theta)
        if (phi, theta) != (0.0, 180.0):
            assert d.T1 == pytest.approx(0.3 / ck, rel=1e-12), (phi, theta)
        # no 1.x rescale: lam_max of a 1.x load would be <= 1; here the 90 deg keys keep T0 = 7.2
    d = ds[(90.0, 90.0)]
    assert lam_max(d.T0, d.T1, d.rho) > 1 + 1e-6


def test_s2_key_without_record(tmp_path):
    """The 2.x record is matched on the raw CSV angles at 1e-9 degree, as BBRHFSSData::LoadWaveguide does: a
    key 1e-7 degree off its declared angle (inside C2's 0.01 degree rounding) has no record, BBR012."""
    def shift(base, stem, phi, theta, dphi):
        for e in (0, 1):
            for name in ("far_field.csv", "waveguide.csv"):
                p = os.path.join(base, f"{stem}_Ephi={e}", name)
                lines = open(p).read().splitlines()
                for i, line in enumerate(lines[1:], start=1):
                    v = line.split(",")
                    if float(v[2]) == phi and float(v[3]) == theta:
                        v[2] = repr(float(v[2]) + dphi)
                        lines[i] = ",".join(v)
                open(p, "w").write("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match=r"^BBR012: .*no per-key record"):
        s2_load(tmp_path / "a", edit=lambda b, st: shift(b, st, 0, 180, 1e-7))     # the dark key: masking not skipped
    with pytest.raises(ValueError, match=r"^BBR012: .*no per-key record"):
        s2_load(tmp_path / "b", edit=lambda b, st: shift(b, st, 45, 135, 1e-7))
    ds = s2_load(tmp_path / "c", edit=lambda b, st: shift(b, st, 45, 135, 1e-10))  # inside the 1e-9 rule
    assert ds[(45.0, 135.0)].T0 > 0


def test_s2_exit_cs_field_frame_refused(tmp_path):
    """CDX-035: the mirror refuses exit-CS field components under 2.x before computing anything, as the C++
    does at sidecar parse; it reads E_t as the HFSS-global (Ex, Ey)."""
    base = str(tmp_path)
    stem = fx.write_nine_key_grid2(base)
    sc = fx.sidecar2()
    sc["exit_field"].update(field_in_ref_cs=True, field_components_frame="exit_cs")
    sidecar.write(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR025: exit_field\.field_components_frame is exit_cs"):
        hfss.load_dataset(stem, base)


def set_rows(key, e, col, value, rows=None):
    """A CSV edit writing `value` into column `col` (4 OutgoingPower, 5 IngoingPower, 15 InterceptFactor) of
    polarization e's rows of one key: every row, or the rows at the indices in `rows`, counted within the key."""
    def edit(base, stem):
        path = os.path.join(base, f"{stem}_Ephi={e}", "waveguide.csv")
        lines, n = open(path).read().splitlines(), 0
        for i, line in enumerate(lines[1:], start=1):
            v = line.split(",")
            if float(v[2]) == key[0] and float(v[3]) == key[1]:
                if rows is None or n in rows:
                    v[col] = value
                    lines[i] = ",".join(v)
                n += 1
        open(path, "w").write("\n".join(lines) + "\n")
    return edit


S2_ROW_EDITS = [   # (id, edit, message): each refused by BBRHFSSData::LoadWaveguide and check_csvs with BBR012
    ("dark-power-first", set_rows((0, 180), 1, 4, "5e-10", [0]), r"\(0, 180\) polarization 1: the sidecar declares this basis dark"),
    ("dark-power-later", set_rows((0, 180), 1, 4, "1e-25", [2]), r"\(0, 180\) polarization 1: the sidecar declares this basis dark"),
    ("negative-first", set_rows((45, 135), 0, 4, "-6e-10", [0]), r"\(45, 135\) polarization 0: OutgoingPower -6e-10 is negative"),
    ("negative-later", set_rows((45, 135), 1, 4, "-1e-12", [3]), r"\(45, 135\) polarization 1: OutgoingPower -1e-12 is negative"),
    ("factor-key", set_rows((45, 135), 0, 15, "0.1"), r"\(45, 135\) polarization 0: InterceptFactor 0\.1 differs from the sidecar"),
    ("factor-later", set_rows((90, 90), 1, 15, "0.0839", [4]), r"\(90, 90\) polarization 1: InterceptFactor 0\.0839 differs from the "
                                                                r"key's first row"),
    ("t-vs-sidecar", set_rows((45, 135), 0, 4, "6.6e-10"), r"\(45, 135\) polarization 0: T = .* but the sidecar's t_flux / c_k"),
    ("ingoing", set_rows((90, 180), 1, 5, "2e-09"), r"\(90, 180\) polarization 1: T = .* but the sidecar's t_flux / c_k"),
]


@pytest.mark.parametrize("edit,message", [x[1:] for x in S2_ROW_EDITS], ids=[x[0] for x in S2_ROW_EDITS])
def test_s2_rows_refused(tmp_path, edit, message):
    """CDX-037: the mirror applies the 2.x row checks of the C++ loader before computing T or masking a
    basis, so it never returns a state the native loader refuses (a dark basis with T > 0, a negative T,
    a T set by a wrong divisor)."""
    with pytest.raises(ValueError, match=r"^BBR012: crack_500GHz_Ephi=\d/waveguide\.csv, key " + message):
        s2_load(tmp_path, edit=edit)


def test_s2_axis_count_one_bbr024(tmp_path):
    """load_dataset on an unchecked sidecar with a one-node rectangle axis: BBR024, as the parse rule says
    (count < 2), not a division by zero in the quadrature."""
    base = str(tmp_path)
    stem = fx.write_nine_key_grid2(base)
    sc = fx.sidecar2()
    sc["exit_field"]["grid"]["y_e"].update({"count": 1, "max": sc["exit_field"]["grid"]["y_e"]["min"], "step": 0.0})
    sidecar.write(base, stem, sc)
    with pytest.raises(ValueError, match=r"^BBR024: .*count 1"):
        hfss.load_dataset(stem, base)
    with pytest.raises(ValueError, match=r"^BBR024"):
        hfss._axis_index(0.0, {"min": 0.0, "max": 0.0, "count": 1}, "Z", (0.0, 180.0))


def test_s2_dark_masking(tmp_path):
    ds = s2_load(tmp_path)
    d = ds[(0.0, 180.0)]
    assert d.dark1 and not d.dark0 and d.T1 == 0.0 and d.rho == 0j
    assert not any((e.dark0 or e.dark1) for k, e in ds.items() if k != (0.0, 180.0))
    assert hfss.transmittance(d, 0.0, 1.0) == 0.0
    assert hfss.transmittance(d, R2, R2) == pytest.approx(0.5 * d.T0, rel=1e-12)   # no cross term
    # the dark basis enters neither sampler: (R2, R2) draws exactly as the bright basis alone
    assert hfss.direction_weights(d, R2, R2) == pytest.approx(hfss.direction_weights(d, 1.0, 0.0), rel=1e-12, abs=0)
    assert hfss.exit_position_weights(d, R2, R2) == pytest.approx(hfss.exit_position_weights(d, 1.0, 0.0), rel=1e-12, abs=0)
    # a tiny bright share still samples the bright basis only, with no BBR027
    assert hfss.direction_weights(d, 1e-9, 1.0) is not None
    assert hfss.exit_position_weights(d, 1e-9, 1.0) == pytest.approx(hfss.exit_position_weights(d, 1.0, 0.0), rel=1e-12, abs=0)
    # the dark basis alone has nothing to sample (the C++ never draws it: T = 0)
    with pytest.raises(ValueError, match="BBR027"):
        hfss.direction_weights(d, 0.0, 1.0)
    with pytest.raises(ValueError, match="BBR027"):
        hfss.exit_position_weights(d, 0.0, 1.0)


def edit_powers_at(key, out0, out1):
    """A CSV edit setting OutgoingPower of every row of one key, per polarization."""
    def edit(base, stem):
        for e, val in ((0, out0), (1, out1)):
            path = os.path.join(base, f"{stem}_Ephi={e}", "waveguide.csv")
            lines = open(path).read().splitlines()
            for i, line in enumerate(lines[1:], start=1):
                v = line.split(",")
                if float(v[2]) == key[0] and float(v[3]) == key[1]:
                    v[4] = repr(val)
                    lines[i] = ",".join(v)
            open(path, "w").write("\n".join(lines) + "\n")
    return edit


def test_s2_rho_weighted(tmp_path):
    # the lower-power variant at (45, 135): the cap is inactive, so rho is visible in the mixed T
    tflux = [[0.6, 0.0 if (ph, t) == (0, 180) else 0.3, 0, 1 if (ph, t) == (0, 180) else 0]
             for ph in (0, 45, 90) for t in (90, 135, 180)]
    tflux[4] = [0.3, 0.15, 0, 0]
    ds = s2_load(tmp_path, tflux=tflux, edit=edit_powers_at((45.0, 135.0), 0.3e-9, 0.15e-9))
    d = ds[(45.0, 135.0)]
    c = np.cos(np.pi / 4)
    T0, T1 = 0.3 / c, 0.15 / c
    assert (d.T0, d.T1) == pytest.approx((T0, T1), rel=1e-12)
    # rho, independently, from the fixture rows: sum_i w_i E_t0,i . conj(E_t1,i) over (Ex, Ey)
    num, p0, p1 = 0j, 0.0, 0.0
    for i in range(9):
        a, b = s2_fields(0, i)[:2], s2_fields(1, i)[:2]
        num += s2_w(i) * np.sum(a * np.conj(b))
        p0 += s2_w(i) * np.sum(np.abs(a) ** 2)
        p1 += s2_w(i) * np.sum(np.abs(b) ** 2)
    rho_w = (num / np.sqrt(p0 * p1)).real
    assert d.rho.real == pytest.approx(rho_w, rel=1e-12, abs=0)
    # the unweighted overlap over (Ex, Ey, Ez), the 1.x rule, is a different number
    n3 = sum(np.sum(s2_fields(0, i) * np.conj(s2_fields(1, i))) for i in range(9))
    q0 = sum(np.sum(np.abs(s2_fields(0, i)) ** 2) for i in range(9))
    q1 = sum(np.sum(np.abs(s2_fields(1, i)) ** 2) for i in range(9))
    rho_u = (n3 / np.sqrt(q0 * q1)).real
    assert abs(rho_w - rho_u) > 1e-4
    mixed = 0.5 * T0 + 0.5 * T1 + np.sqrt(T0 * T1) * rho_w
    assert 0.0 < mixed < 1.0                                    # the cap is inactive
    got = hfss.transmittance(d, R2, R2)
    assert got == pytest.approx(mixed, abs=1e-12)
    assert abs(got - (0.5 * T0 + 0.5 * T1 + np.sqrt(T0 * T1) * rho_u)) > 1e-5
    assert got == pytest.approx(0.6050564495, abs=1e-9)         # the number the C++ test sees


def test_s2_mixture_clamped(tmp_path):
    ds = s2_load(tmp_path)
    d = ds[(0.0, 90.0)]
    psi = (np.arange(360) + 0.5) / 360 * 2 * np.pi
    raw = np.cos(psi) ** 2 * d.T0 + np.sin(psi) ** 2 * d.T1 + 2 * np.cos(psi) * np.sin(psi) * np.sqrt(d.T0 * d.T1) * d.rho.real
    clamped = np.clip(raw, 0.0, 1.0)
    mean, w = hfss.random_polarization_mixture(d)
    assert mean == pytest.approx(clamped.mean(), abs=1e-12)
    assert mean < (d.T0 + d.T1) / 2 and mean <= 1.0
    assert w.sum() == pytest.approx(1.0)


def test_schema2_header_without_sidecar_refused(tmp_path):
    base = str(tmp_path)
    stem = fx.write_nine_key_grid2(base)       # 16-column waveguide.csv, no sidecar
    with pytest.raises(ValueError, match="^BBR024"):
        hfss.load_dataset(stem, base)


def test_schema2_sidecar_with_1x_header_refused(tmp_path):
    base = str(tmp_path / "s1")
    write_dataset(tmp_path / "s1", "c_500GHz")  # a 1.0-header dataset
    sc = fx.sidecar2(stem="c_500GHz")
    sidecar.write(base, "c_500GHz", sc)
    with pytest.raises(ValueError, match="BBR013"):
        hfss.load_dataset("c_500GHz", base)


def test_s1_mirror_unchanged():
    base = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "waveguides")
    pinned = {"InfParallelPlate_crack1Rohan_500GHz": (0.5143466820199462, 0.4858421150130238,
                                                        0.9996221702766228 + 0.001069073061348357j),
              "InfParallelPlate_crack2_500GHz": (0.5015062206178313, 0.4985601485973919,
                                                   0.9998672692267748 - 6.399260932677592e-05j)}
    for stem, (T0, T1, rho) in pinned.items():
        ds = hfss.load_dataset(stem, base)
        d = ds[(45.0, 180.0)]
        assert d.T0 == pytest.approx(T0, abs=1e-15, rel=0) and d.T1 == pytest.approx(T1, abs=1e-15, rel=0)
        assert abs(d.rho - rho) <= 1e-15
        assert d.schema == 1 and not d.dark0 and not d.dark1 and np.all(d.w == 1.0)
