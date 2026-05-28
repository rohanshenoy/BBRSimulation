#!/usr/bin/env python3
"""
plot_diffraction.py — compare Geant4 diffraction output vs HFSS reference.
Run from the build/ directory after: ./BBRSim diffraction.mac

Crack-local frame (standard geometry, no rotation):
  normal_hat = x_world  (exit face outward normal / propagation axis)
  theta_hat  = y_world  (long dimension)
  phi_hat    = z_world  (gap dimension b)

Outgoing direction formula (from SampleOutgoingDirection):
  dir_out = sin(T)*cos(P)*normal_hat + sin(T)*sin(P)*theta_hat + cos(T)*phi_hat
so: T = acos(dir_z),  P = atan2(dir_y, dir_x)

Exit position (from SampleExitPosition, crack centred at origin):
  pos_out.y / m  ≈  HFSS Y  (long dimension, metres)
  pos_out.z / m  ≈  HFSS Z  (gap dimension, metres)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA = '../HFSSSimData/InfParallelPlate_crack1Rohan_500GHz'

# ── Geant4 output ─────────────────────────────────────────────────────────────
try:
    g4 = pd.read_csv('diffraction_output.csv')
except FileNotFoundError:
    sys.exit('diffraction_output.csv not found — run ./BBRSim diffraction.mac first')

g4['theta'] = np.degrees(np.arccos(np.clip(g4['dir_z'], -1.0, 1.0)))
g4['phi']   = np.degrees(np.arctan2(g4['dir_y'], g4['dir_x']))
print(f'Transmitted photons recorded: {len(g4)}')

# ── HFSS CSVs — normal incidence filter ──────────────────────────────────────
# Normal incidence → IWaveTheta = 180° (photon k-vector anti-parallel to ẑ_i=normal_hat,
# i.e. propagating from inside crack toward exit face). IWavePhi = 0°.
def load_and_filter(ephi):
    d = f'{DATA}_Ephi={ephi}'
    ff = pd.read_csv(f'{d}/far_field.csv')
    wg = pd.read_csv(f'{d}/waveguide.csv')
    mask_ff = (ff['IWavePhi'].abs() < 0.5) & (ff['IWaveTheta'] > 179.5)
    mask_wg = (wg['IWavePhi'].abs() < 0.5) & (wg['IWaveTheta'] > 179.5)
    return ff[mask_ff].reset_index(drop=True), wg[mask_wg].reset_index(drop=True)

ff0, wg0 = load_and_filter(0)
ff1, wg1 = load_and_filter(1)
print(f'HFSS far-field rows (normal inc): {len(ff0)},  waveguide rows: {len(wg0)}')

# 45/45 polarization: E_theta = E_phi = 1/√2
Et = Ep = 1.0 / np.sqrt(2.0)

# ── Far-field combined power for 45/45 pol ───────────────────────────────────
Eth_re = Et*ff0['rEtheta_real'].values + Ep*ff1['rEtheta_real'].values
Eth_im = Et*ff0['rEtheta_imag'].values + Ep*ff1['rEtheta_imag'].values
Eph_re = Et*ff0['rEphi_real'].values   + Ep*ff1['rEphi_real'].values
Eph_im = Et*ff0['rEphi_imag'].values   + Ep*ff1['rEphi_imag'].values
ff0['power'] = Eth_re**2 + Eth_im**2 + Eph_re**2 + Eph_im**2
ff0['Theta'] = ff0['Theta']  # already present
ff0['Phi']   = ff0['Phi']

# ── Waveguide combined |E|² for 45/45 pol ────────────────────────────────────
Px_re = Et*wg0['Ex_real'].values + Ep*wg1['Ex_real'].values
Px_im = Et*wg0['Ex_imag'].values + Ep*wg1['Ex_imag'].values
Py_re = Et*wg0['Ey_real'].values + Ep*wg1['Ey_real'].values
Py_im = Et*wg0['Ey_imag'].values + Ep*wg1['Ey_imag'].values
Pz_re = Et*wg0['Ez_real'].values + Ep*wg1['Ez_real'].values
Pz_im = Et*wg0['Ez_imag'].values + Ep*wg1['Ez_imag'].values
wg0 = wg0.copy()
wg0['power'] = Px_re**2+Px_im**2+Py_re**2+Py_im**2+Pz_re**2+Pz_im**2

# Sum over Z (gap dimension) → 1-D distribution in Y (long dimension)
wg_y = wg0.groupby('Y')['power'].sum().reset_index().sort_values('Y')

# ── Helper: normalise HFSS curve to match Geant4 histogram area ──────────────
def norm_to_hist(x_hfss, p_hfss, bin_edges, n_g4):
    """Return HFSS (x, y) scaled so its integral matches the G4 histogram."""
    bin_w   = bin_edges[1] - bin_edges[0]
    centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    p_interp = np.interp(centers, x_hfss, p_hfss, left=0, right=0)
    scale = (n_g4.sum() * bin_w) / (p_interp.sum() * bin_w) if p_interp.sum() > 0 else 1
    return x_hfss, p_hfss * scale

# ── Plot 1: outgoing phi at Theta ≈ 90° ──────────────────────────────────────
theta_tol = 1.5   # degrees — matches how densely sampled the CSV is
g4_t90  = g4[np.abs(g4['theta'] - 90.0) < theta_tol]
ff_t90  = ff0[np.abs(ff0['Theta'] - 90.0) < theta_tol].sort_values('Phi')
print(f'G4 events with Theta in [89,91] deg: {len(g4_t90)}')

bins1   = np.linspace(-30, 30, 41)
n1, _   = np.histogram(g4_t90['phi'], bins=bins1)

fig1, ax1 = plt.subplots(figsize=(8, 5))
ax1.hist(g4_t90['phi'], bins=bins1, label='Geant4', color='steelblue', alpha=0.7)

if len(ff_t90) > 1:
    xh, yh = norm_to_hist(ff_t90['Phi'].values, ff_t90['power'].values, bins1, n1)
    ax1.plot(xh, yh, 'r-', lw=2, label='HFSS')

ax1.set_xlabel('Outgoing phi (degrees)')
ax1.set_ylabel('Counts')
ax1.set_title(f'Output phi Distribution at Theta=90 degrees, number of events = {len(g4_t90)}')
ax1.legend()
plt.tight_layout()
plt.savefig('plot_phi_dist.png', dpi=150)
print('Saved plot_phi_dist.png')

# ── Plot 2: exit-face Y position (long dimension) ────────────────────────────
bins2   = np.linspace(-0.006, 0.006, 61)
n2, _   = np.histogram(g4['pos_y_m'], bins=bins2)

fig2, ax2 = plt.subplots(figsize=(8, 5))
ax2.hist(g4['pos_y_m'], bins=bins2, label='Geant4', color='sandybrown', alpha=0.7)

xh2, yh2 = norm_to_hist(wg_y['Y'].values, wg_y['power'].values, bins2, n2)
ax2.plot(xh2, yh2, 'b-', lw=2, label='HFSS')

ax2.set_xlabel('Output Y Position (m)  [long dimension]')
ax2.set_ylabel('Counts')
ax2.set_title('Exit Face Position Distribution (long dimension)')
ax2.legend()
plt.tight_layout()
plt.savefig('plot_pos_dist.png', dpi=150)
print('Saved plot_pos_dist.png')

plt.show()
