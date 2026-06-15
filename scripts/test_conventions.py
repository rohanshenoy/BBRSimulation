#!/usr/bin/env python3
"""
test_conventions.py — unit tests for BBRSim HFSS angle / basis conventions.

Layer 1: Pure math, no CSV files needed.
Layer 2: HFSS data contract, needs HFSSSimData/ relative to this file.

Run from repo root:
  python3 test_conventions.py
"""

import sys, math, pathlib
import numpy as np

# ---------------------------------------------------------------------------
# Crack-local basis (aligned-crack geometry only — these are hardcoded in
# BBSimOpBoundaryProcess.cc for the current test geometry).
# ---------------------------------------------------------------------------
normal_hat = np.array([1., 0., 0.])   # +x_world  exit face outward normal
theta_hat  = np.array([0., 1., 0.])   # +y_world  long dimension
phi_hat    = np.array([0., 0., 1.])   # +z_world  gap dimension b

# ---------------------------------------------------------------------------
# Replicate the C++ formulas exactly (BBSimOpBoundaryProcess + BBRHFSSData)
# ---------------------------------------------------------------------------

def iwave_angles(khat):
    """Return (IWaveTheta_deg, IWavePhi_raw_deg, IWavePhi_folded_deg)."""
    khat = np.asarray(khat, float)
    cosVal = float(np.clip(-np.dot(khat, normal_hat), -1., 1.))
    theta  = math.degrees(math.acos(cosVal))
    phi_raw = math.degrees(math.atan2(-np.dot(khat, theta_hat),
                                       np.dot(khat, phi_hat)))
    phi_fold = abs(phi_raw)
    if phi_fold > 90.:
        phi_fold = 180. - phi_fold
    return theta, phi_raw, phi_fold


def khat_from_iwave(T_deg, P_raw_deg):
    """Inverse of iwave_angles: propagation direction khat from (T, P_raw).

    Derivation: IWaveTheta = acos(-khat·z_i) where z_i = normal_hat.
    IWavePhi = atan2(khat·y_i, khat·x_i) where x_i=phi_hat, y_i=-theta_hat.
    So khat·z_i = -cos(T), khat·x_i = sin(T)*cos(P), khat·y_i = sin(T)*sin(P).
    """
    T = math.radians(T_deg)
    P = math.radians(P_raw_deg)
    return (-math.cos(T) * normal_hat
            - math.sin(T) * math.sin(P) * theta_hat
            + math.sin(T) * math.cos(P) * phi_hat)


def incoming_basis(T_deg, P_raw_deg):
    """eTheta_in, ePhi_in — exactly as in HandleDiffractionBoundary."""
    th = math.radians(T_deg)
    ph = math.radians(P_raw_deg)
    eT = (+ math.sin(th) * normal_hat
          - math.cos(th) * math.sin(ph) * theta_hat
          + math.cos(th) * math.cos(ph) * phi_hat)
    eP = (- math.cos(ph) * theta_hat
          - math.sin(ph) * phi_hat)
    return eT, eP


def outgoing_dir(T_deg, P_deg):
    """dir = sinT*cosP*normal + sinT*sinP*theta + cosT*phi (SampleOutgoingDirection)."""
    T = math.radians(T_deg)
    P = math.radians(P_deg)
    return (math.sin(T)*math.cos(P) * normal_hat
          + math.sin(T)*math.sin(P) * theta_hat
          + math.cos(T)             * phi_hat)


def outgoing_basis(T_deg, P_deg):
    """eTh_out, ePhi_out — spherical basis for outgoing frame."""
    T = math.radians(T_deg)
    P = math.radians(P_deg)
    eTh = (  math.cos(T)*math.cos(P) * normal_hat
           + math.cos(T)*math.sin(P) * theta_hat
           - math.sin(T)             * phi_hat)
    ePh = (- math.sin(P) * normal_hat
           + math.cos(P) * theta_hat)
    return eTh, ePh


# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------
PASS = FAIL = 0

def check(name, cond, detail=''):
    global PASS, FAIL
    if cond:
        print(f'  PASS  {name}')
        PASS += 1
    else:
        print(f'  FAIL  {name}  [{detail}]')
        FAIL += 1

def near(a, b, tol=1e-9):
    return abs(float(a) - float(b)) < tol

def vec_near(a, b, tol=1e-9):
    return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float))) < tol

# ---------------------------------------------------------------------------
# Layer 1a: IWaveTheta / IWavePhi
# ---------------------------------------------------------------------------
print('=== Layer 1a: IWaveTheta ===')

T, _, _ = iwave_angles([1, 0, 0])
check('khat=+normal_hat  → IWaveTheta=180°', near(T, 180.),
      f'got {T:.4f}')

T, _, _ = iwave_angles([-1, 0, 0])
check('khat=-normal_hat  → IWaveTheta=0°',   near(T, 0.),
      f'got {T:.4f}')

T, _, _ = iwave_angles([0, 0, 1])
check('khat=+phi_hat     → IWaveTheta=90°',  near(T, 90.),
      f'got {T:.4f}')

T, _, _ = iwave_angles([0, 1, 0])
check('khat=+theta_hat   → IWaveTheta=90°',  near(T, 90.),
      f'got {T:.4f}')

T, _, _ = iwave_angles([1/math.sqrt(2), 0, 1/math.sqrt(2)])
check('khat 45° in x-z  → IWaveTheta=135°', near(T, 135., 1e-6),
      f'got {T:.4f}')

print()
print('=== Layer 1a: IWavePhi ===')

_, P_raw, P_fold = iwave_angles([1, 0, 0])
check('khat=+normal_hat  → IWavePhi=0°',  near(P_fold, 0.),
      f'got {P_fold:.4f}')

_, P_raw, P_fold = iwave_angles([0, 0, 1])
check('khat=+phi_hat     → IWavePhi=0°',  near(P_fold, 0.),
      f'got {P_fold:.4f}  (phi_hat is x_i axis, Phi=0°)')

_, P_raw, P_fold = iwave_angles([0, 1, 0])
check('khat=+theta_hat   → IWavePhi=90°', near(P_fold, 90.),
      f'got {P_fold:.4f}')

_, P_raw, P_fold = iwave_angles([0, -1, 0])
check('khat=-theta_hat   → IWavePhi=90°  (quarter-sym fold)', near(P_fold, 90.),
      f'got {P_fold:.4f}')

# ---------------------------------------------------------------------------
# Layer 1b: Incoming basis — orthonormality and perpendicularity to khat
# ---------------------------------------------------------------------------
print()
print('=== Layer 1b: Incoming basis (eTheta_in, ePhi_in) ===')

test_iwave = [
    (180., 0.,  'normal incidence'),
    (135., 0.,  '45° in x-z plane'),
    (135., 45., '45° off-axis'),
    ( 90., 0.,  'grazing, P=0'),
    ( 90., 45., 'grazing, P=45'),
    ( 45., 30., 'oblique'),
]

for T_deg, P_raw, label in test_iwave:
    eT, eP = incoming_basis(T_deg, P_raw)
    khat   = khat_from_iwave(T_deg, P_raw)

    check(f'{label:25s}  |eTheta|=1',     near(np.linalg.norm(eT), 1., 1e-9),
          f'{np.linalg.norm(eT):.2e}')
    check(f'{label:25s}  |ePhi|=1',       near(np.linalg.norm(eP), 1., 1e-9),
          f'{np.linalg.norm(eP):.2e}')
    check(f'{label:25s}  eTheta⊥ePhi',    near(np.dot(eT, eP), 0., 1e-9),
          f'{np.dot(eT, eP):.2e}')
    check(f'{label:25s}  eTheta⊥khat',    near(np.dot(eT, khat), 0., 1e-9),
          f'{np.dot(eT, khat):.2e}')
    check(f'{label:25s}  ePhi⊥khat',      near(np.dot(eP, khat), 0., 1e-9),
          f'{np.dot(eP, khat):.2e}')

# ---------------------------------------------------------------------------
# Layer 1c: Polarization decomposition — test case (45/45 pol)
# ---------------------------------------------------------------------------
print()
print('=== Layer 1c: Polarization decomposition (diffraction.mac test case) ===')

# phat = (0, -1/√2, -1/√2) → should give E_theta = E_phi = 1/√2
phat_test = np.array([0., -1./math.sqrt(2.), -1./math.sqrt(2.)])
T_ni, P_raw_ni, _ = iwave_angles([1, 0, 0])   # normal incidence
eT_ni, eP_ni = incoming_basis(T_ni, P_raw_ni)

E_theta = np.dot(phat_test, eT_ni)
E_phi   = np.dot(phat_test, eP_ni)
inv_sqrt2 = 1. / math.sqrt(2.)

check('normal inc, 45/45 pol:  E_theta = 1/√2', near(E_theta, inv_sqrt2, 1e-9),
      f'got {E_theta:.6f}')
check('normal inc, 45/45 pol:  E_phi   = 1/√2', near(E_phi,   inv_sqrt2, 1e-9),
      f'got {E_phi:.6f}')
check('normal inc, 45/45 pol:  |phat|=1 consistency',
      near(E_theta**2 + E_phi**2, 1., 1e-9),
      f'E²={E_theta**2 + E_phi**2:.6f}')

# ---------------------------------------------------------------------------
# Layer 1d: Outgoing direction formula — anchor points
# ---------------------------------------------------------------------------
print()
print('=== Layer 1d: Outgoing direction (T,P → world dir) ===')

check('T=90, P=0    → +normal_hat (forward)',  vec_near(outgoing_dir(90,  0), normal_hat))
check('T=90, P=90   → +theta_hat  (long dim)', vec_near(outgoing_dir(90, 90), theta_hat))
check('T=90, P=-90  → -theta_hat',             vec_near(outgoing_dir(90,-90),-theta_hat))
check('T=0          → +phi_hat    (gap end-on)',vec_near(outgoing_dir( 0,  0), phi_hat))
check('T=180        → -phi_hat',               vec_near(outgoing_dir(180, 0),-phi_hat))

# ---------------------------------------------------------------------------
# Layer 1e: Outgoing basis — orthonormality and perpendicularity to dir
# ---------------------------------------------------------------------------
print()
print('=== Layer 1e: Outgoing basis (eTh_out, ePhi_out) ===')

test_out = [(90, 0), (90, 45), (45, 30), (60, -60), (30, 90)]
for T_deg, P_deg in test_out:
    d   = outgoing_dir(T_deg, P_deg)
    eTh, ePh = outgoing_basis(T_deg, P_deg)

    check(f'T={T_deg:3d},P={P_deg:4d}  |eTh|=1',    near(np.linalg.norm(eTh), 1., 1e-9))
    check(f'T={T_deg:3d},P={P_deg:4d}  |ePh|=1',    near(np.linalg.norm(ePh), 1., 1e-9))
    check(f'T={T_deg:3d},P={P_deg:4d}  eTh⊥ePh',   near(np.dot(eTh, ePh), 0., 1e-9))
    check(f'T={T_deg:3d},P={P_deg:4d}  eTh⊥dir',   near(np.dot(eTh, d),   0., 1e-9))
    check(f'T={T_deg:3d},P={P_deg:4d}  ePh⊥dir',   near(np.dot(ePh, d),   0., 1e-9))

# ---------------------------------------------------------------------------
# Layer 2: HFSS data contract
# ---------------------------------------------------------------------------
print()
print('=== Layer 2: HFSS data contract ===')

DATA = pathlib.Path(__file__).parent / 'HFSSSimData' / 'InfParallelPlate_crack1Rohan_500GHz'

try:
    import pandas as pd
except ImportError:
    print('  SKIP  (pandas not available)')
    pd = None

if pd is not None:
    try:
        wg0 = pd.read_csv(str(DATA) + '_Ephi=0/waveguide.csv')
        wg1 = pd.read_csv(str(DATA) + '_Ephi=1/waveguide.csv')
        ff0 = pd.read_csv(str(DATA) + '_Ephi=0/far_field.csv')
        ff1 = pd.read_csv(str(DATA) + '_Ephi=1/far_field.csv')

        # Transmittance at normal incidence (IWaveTheta=180, IWavePhi=0)
        def T_at(wg, iwT, iwP_tol=0.5):
            sub = wg[(wg['IWaveTheta'].between(iwT-0.5, iwT+0.5)) &
                     (wg['IWavePhi'].abs() < iwP_tol)]
            row = sub.iloc[0]
            return row['OutgoingPower'] / row['IngoingPower']

        T_ni_e0 = T_at(wg0, 180.)
        T_ni_e1 = T_at(wg1, 180.)
        T_zero_e0 = T_at(wg0, 0.)

        check('Ephi=0, IWaveTheta=180: T > 1.0   (confirmed >1)',
              T_ni_e0 > 1.0, f'T={T_ni_e0:.4f}')
        check('Ephi=1, IWaveTheta=180: T ≈ 0',
              T_ni_e1 < 1e-6, f'T={T_ni_e1:.2e}')
        check('Ephi=0, IWaveTheta=0:   T ≈ 0     (backward, blocked)',
              T_zero_e0 < 1e-20, f'T={T_zero_e0:.2e}')

        # 45/45 pol transmittance: E²·T0 + E²·T1
        Et = Ep = inv_sqrt2
        T_combined = Et**2 * T_ni_e0 + Ep**2 * T_ni_e1
        check('45/45 pol: combined T ≈ 0.528 (± 0.05)',
              abs(T_combined - 0.528) < 0.05, f'T={T_combined:.4f}')

        # Far-field: dominant power at T=90 (forward direction, after sinT weighting)
        sub0 = ff0[(ff0['IWaveTheta'] > 179.5) & (ff0['IWavePhi'].abs() < 0.5)].reset_index(drop=True)
        sub1 = ff1[(ff1['IWaveTheta'] > 179.5) & (ff1['IWavePhi'].abs() < 0.5)].reset_index(drop=True)

        Eth_re = Et*sub0['rEtheta_real'].values + Ep*sub1['rEtheta_real'].values
        Eth_im = Et*sub0['rEtheta_imag'].values + Ep*sub1['rEtheta_imag'].values
        Eph_re = Et*sub0['rEphi_real'].values   + Ep*sub1['rEphi_real'].values
        Eph_im = Et*sub0['rEphi_imag'].values   + Ep*sub1['rEphi_imag'].values
        power  = Eth_re**2 + Eth_im**2 + Eph_re**2 + Eph_im**2

        T_rad = np.radians(sub0['Theta'].values)
        P_rad = np.radians(sub0['Phi'].values)
        sinT  = np.sin(T_rad)
        w     = power * sinT            # solid-angle weighted
        dir_x = np.sin(T_rad)*np.cos(P_rad)   # along normal_hat
        mean_x = (w * dir_x).sum() / w.sum()

        check('sinT-weighted mean dir_x > 0.6 (mostly forward)',
              mean_x > 0.6, f'mean_x={mean_x:.3f}')

        # Waveguide exit positions are within crack dimensions
        sub_wg = wg0[(wg0['IWaveTheta'] > 179.5) & (wg0['IWavePhi'].abs() < 0.5)]
        half_b = 25e-6   # b/2 = 25 µm
        half_a = 5e-3    # a/2 = 5 mm (long dim)
        check('Waveguide Z within ±b/2 = ±25 µm',
              sub_wg['Z'].abs().max() <= half_b + 1e-9,
              f'max |Z|={sub_wg["Z"].abs().max():.2e} m')
        check('Waveguide Y within ±a/2 = ±5 mm',
              sub_wg['Y'].abs().max() <= half_a + 1e-9,
              f'max |Y|={sub_wg["Y"].abs().max():.2e} m')

    except FileNotFoundError as e:
        print(f'  SKIP  CSV not found: {e}')

# ---------------------------------------------------------------------------
# Layer 3: Exit-position CDF — dimension mapping
# ---------------------------------------------------------------------------
print()
print('=== Layer 3: Exit-position CDF dimension mapping ===')

if pd is not None:
    try:
        HALF_A = 5e-3    # m  long dim a = 10 mm
        HALF_B = 25e-6   # m  gap dim b = 50 µm

        ni0 = wg0[(wg0['IWaveTheta'] > 179.5) & (wg0['IWavePhi'].abs() < 0.5)].reset_index(drop=True)
        ni1 = wg1[(wg1['IWaveTheta'] > 179.5) & (wg1['IWavePhi'].abs() < 0.5)].reset_index(drop=True)

        # --- 3a: CSV column identity — which physical dimension is Y vs Z ---
        check('CSV X = 0 everywhere (exit face, not propagation interior)',
              (ni0['X'] == 0.0).all(),
              f'X unique: {ni0["X"].unique()}')
        check('CSV Y spans ±a/2 = ±5 mm  (long dimension)',
              abs(ni0['Y'].max() - HALF_A) < 1e-5 and abs(ni0['Y'].min() + HALF_A) < 1e-5,
              f'Y=[{ni0["Y"].min():.4e}, {ni0["Y"].max():.4e}] m')
        check('CSV Z spans ±b/2 = ±25 µm (gap dimension)',
              abs(ni0['Z'].max() - HALF_B) < 1e-7 and abs(ni0['Z'].min() + HALF_B) < 1e-7,
              f'Z=[{ni0["Z"].min():.2e}, {ni0["Z"].max():.2e}] m')
        y_span = ni0['Y'].max() - ni0['Y'].min()
        z_span = ni0['Z'].max() - ni0['Z'].min()
        check('Y span / Z span > 100  (long dim >> gap, column roles unambiguous)',
              y_span / z_span > 100,
              f'ratio={y_span / z_span:.1f}')

        # --- 3b: Y-power distribution shape ---
        # Build combined |E|² for 45/45 polarization at normal incidence.
        Et45 = Ep45 = 1.0 / math.sqrt(2.0)
        ni0 = ni0.copy()
        ni0['pw'] = (
            (Et45*ni0['Ex_real'] + Ep45*ni1['Ex_real'])**2 +
            (Et45*ni0['Ex_imag'] + Ep45*ni1['Ex_imag'])**2 +
            (Et45*ni0['Ey_real'] + Ep45*ni1['Ey_real'])**2 +
            (Et45*ni0['Ey_imag'] + Ep45*ni1['Ey_imag'])**2 +
            (Et45*ni0['Ez_real'] + Ep45*ni1['Ez_real'])**2 +
            (Et45*ni0['Ez_imag'] + Ep45*ni1['Ez_imag'])**2
        )
        wy = ni0.groupby('Y')['pw'].sum()
        peak_Y = abs(wy.index[wy.values.argmax()])
        check('Y-power peak is at |Y| > 3 mm (edge-peaked, not centre-peaked)',
              peak_Y > 3e-3, f'peak |Y|={peak_Y*1e3:.2f} mm')
        # Power-weighted mean |Y| from CSV (used below to validate sampling).
        mean_abs_Y_theory = (np.abs(wy.index.values) * wy.values).sum() / wy.values.sum()

        # --- 3c: Sample from CDF; verify spread in ep_y vs ep_z ---
        # A Y↔Z axis swap would give std_y ≈ 10 µm and std_z ≈ 3 mm — catches
        # both swapped CSV columns and swapped crack_x/crack_y call arguments.
        rng = np.random.default_rng(42)
        N_SAMP = 5000
        cdf_v = np.cumsum(ni0['pw'].values)
        cdf_v /= cdf_v[-1]
        U = rng.uniform(0, 1, N_SAMP)
        idx = np.searchsorted(cdf_v, U).clip(0, len(cdf_v) - 1)
        ep_y = ni0['Y'].values[idx]   # m, CSV Y → long dimension
        ep_z = ni0['Z'].values[idx]   # m, CSV Z → gap dimension

        std_y = ep_y.std()
        std_z = ep_z.std()
        check('Sampled ep_y std > 2 mm  (matches long-dim spread)',
              std_y > 2e-3,  f'std_y={std_y*1e3:.2f} mm')
        check('Sampled ep_z std < 15 µm (matches gap spread)',
              std_z < 15e-6, f'std_z={std_z*1e6:.2f} µm')
        check('ep_y/ep_z std ratio > 100 (Y and Z axes not swapped in CDF)',
              std_y / std_z > 100, f'ratio={std_y / std_z:.1f}')
        check('Sampled ep_y within ±a/2 = ±5 mm',
              ep_y.max() <= HALF_A + 1e-6 and ep_y.min() >= -HALF_A - 1e-6,
              f'ep_y=[{ep_y.min()*1e3:.2f}, {ep_y.max()*1e3:.2f}] mm')
        check('Sampled ep_z within ±b/2 = ±25 µm',
              ep_z.max() <= HALF_B + 1e-7 and ep_z.min() >= -HALF_B - 1e-7,
              f'ep_z=[{ep_z.min()*1e6:.2f}, {ep_z.max()*1e6:.2f}] µm')
        # Sampled mean |ep_y| should match the power-weighted expectation from the CSV
        # (within ~0.3 mm statistical tolerance for 5000 samples).
        mean_abs_Y_sampled = abs(ep_y).mean()
        check('Sampled mean |ep_y| matches CSV power-weighted mean (within 0.3 mm)',
              abs(mean_abs_Y_sampled - mean_abs_Y_theory) < 3e-4,
              f'sampled={mean_abs_Y_sampled*1e3:.2f} mm, theory={mean_abs_Y_theory*1e3:.2f} mm')

        # --- 3d: World-frame mapping (replicate SampleExitPosition) ---
        # C++ code: pos_out = center + (ep.y * CLHEP::m)*crack_x + (ep.z * CLHEP::m)*crack_y
        #           crack_x = theta_hat = [0,1,0], crack_y = phi_hat = [0,0,1]
        # CLHEP::m = 1000 (Geant4 base unit is mm; CSV is in metres).
        # If crack_x/crack_y were swapped the world-y range would shrink to ±0.025 mm.
        # If CLHEP::m were missing the world-y range would shrink to ±0.005 mm.
        CLHEP_m = 1000.0  # mm per metre
        pos_y = ep_y * CLHEP_m   # mm, displacement along theta_hat (world y)
        pos_z = ep_z * CLHEP_m   # mm, displacement along phi_hat  (world z)
        check('World pos_y range > ±3 mm (long dim → theta_hat mapping OK)',
              pos_y.max() > 3.0 and pos_y.min() < -3.0,
              f'pos_y=[{pos_y.min():.2f}, {pos_y.max():.2f}] mm')
        check('World |pos_z| < 0.030 mm (gap → phi_hat mapping OK)',
              abs(pos_z).max() < 0.030,
              f'max|pos_z|={abs(pos_z).max()*1e3:.3f} µm')
        check('World pos_y/pos_z std ratio > 100 (not swapped in world frame)',
              pos_y.std() / pos_z.std() > 100,
              f'ratio={pos_y.std() / pos_z.std():.1f}')

    except FileNotFoundError as e:
        print(f'  SKIP  Layer 3: {e}')
    except Exception as e:
        import traceback
        print(f'  ERROR  Layer 3: {e}')
        traceback.print_exc()

# ---------------------------------------------------------------------------
print()
print(f'Results: {PASS} passed, {FAIL} failed')
sys.exit(0 if FAIL == 0 else 1)
