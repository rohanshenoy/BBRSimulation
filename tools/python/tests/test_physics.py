"""bbrsim.physics: pure functions, no I/O."""
import warnings

import numpy as np
import pytest

from bbrsim import hfss, physics


# --- rows moved from validation/check_physics.py -------------------------------
# Characterization, not a reference: 4.9e-5 / 1.0e-3 / 6.3e-4 (the D values the
# user guide quotes, to two figures) are outputs of this same Drude model, not
# measurements. The 2 % band checks that the quoted numbers still describe it.
@pytest.mark.parametrize("rrr,expected", [(100, 4.9e-5), (3, 1.0e-3), (6, 6.3e-4)])
def test_drude_documented_values_500GHz_4K(rrr, expected):
    assert physics.drude_absorptance(500e9, rrr, 4.0) == pytest.approx(expected, rel=0.02)


def test_hagen_rubens_vs_serov_of_cu():
    assert physics.hagen_rubens_absorptance(150e9, 1.0 / 0.56e-8) == pytest.approx(0.58e-3, rel=0.10)


def test_planck_peak():
    T = 4.0
    u = np.linspace(0.05, 8.0, 4000)
    pdf = physics.planck_photon_number_pdf(u * physics.K_EV * T, T)
    assert u[np.argmax(pdf)] == pytest.approx(physics.PLANCK_PEAK_U, rel=0.02)


# --- new properties no validator states ---------------------------------------
def test_drude_tends_to_hagen_rubens_when_omega_tau_small():
    # RRR 1 at 4 K: tau ~ 2.5e-14 s, f_break ~ 6 THz; 1 GHz is deep in omega*tau << 1.
    f = 1e9
    d = physics.drude_absorptance(f, 1, 4.0)
    hr = physics.hagen_rubens_absorptance(f, physics.sigma_dc(1, 4.0))
    assert d == pytest.approx(hr, rel=0.02)


def test_relaxation_plateau_rrr100():
    # D -> 2/(omega_p tau) above f_break (omega_p = 1.64e16 rad/s for Cu).
    tau = physics.drude_tau(100, 4.0)
    plateau = 2.0 / (1.64e16 * tau)
    assert physics.drude_absorptance(3e12, 100, 4.0) == pytest.approx(plateau, rel=0.05)


def test_hagen_rubens_overestimates_in_relaxation_regime():
    f = 500e9
    assert physics.hagen_rubens_absorptance(f, physics.sigma_dc(100)) > 2 * physics.drude_absorptance(f, 100)


def test_reflectance_bounded_and_vectorised():
    f = np.logspace(10, np.log10(20e12), 24)
    R = physics.drude_reflectance(f, 100, 4.0)
    assert R.shape == f.shape and np.all((R >= 0) & (R <= 1))


def test_absorptance_decreases_with_rrr():
    d = [physics.drude_absorptance(500e9, r, 4.0) for r in (1, 3, 6, 10, 100, 500)]
    assert all(a > b for a, b in zip(d, d[1:]))


def test_sigma_dc_branches():
    assert physics.sigma_dc(100, 4.0) == 100 * physics.SIGMA_RT
    assert physics.sigma_dc(100, 49.999) == 100 * physics.SIGMA_RT
    s50 = physics.sigma_dc(100, 50.0)
    assert s50 == pytest.approx(1 / (1 / (100 * physics.SIGMA_RT) + 1 / (physics.SIGMA_RT * 273 / 50)))


def test_sigma_dc_vectorised_in_T():
    T = np.array([4.0, 49.999, 50.0, 300.0])
    s = physics.sigma_dc(100, T)
    assert isinstance(s, np.ndarray) and s.shape == T.shape
    assert s.tolist() == [physics.sigma_dc(100, float(t)) for t in T]
    assert isinstance(physics.sigma_dc(100, 4.0), float)


def test_sigma_dc_step_at_50K_is_documented_behaviour():
    # Pins the current step: RRR 100 conductivity drops ~19x across 50 K.
    ratio = physics.sigma_dc(100, 49.999) / physics.sigma_dc(100, 50.0)
    assert ratio == pytest.approx(19.3, rel=0.01)


def test_planck_pdf_edges_quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        v = physics.planck_photon_number_pdf(np.array([-1.0, 0.0, 1e-6, 10.0]), 4.0)
    assert v[0] == 0 and v[1] == 0 and v[2] > 0 and v[3] == 0


def test_ev_hz_round_trip():
    f = np.array([10e9, 500e9, 20e12])
    assert np.allclose(physics.ev_to_hz(physics.hz_to_ev(f)), f, rtol=1e-15)


def test_planck_weighted_absorptance_within_range():
    f = np.logspace(10, np.log10(20e12), 400)
    d = physics.drude_absorptance(f, 100, 4.0)
    A = physics.planck_weighted_absorptance(100, 4.0)
    assert d.min() <= A <= d.max()


def test_planck_constants_agree_between_modules():
    # physics.H_EVS mirrors BBRMaterials.hh h_eVs; hfss.H_EV_S mirrors CLHEP.
    assert physics.H_EVS == pytest.approx(hfss.H_EV_S, rel=1e-9)
