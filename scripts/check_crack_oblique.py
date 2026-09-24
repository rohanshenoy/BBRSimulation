"""
check_crack_oblique.py
Validate the HFSS azimuth fold/unfold at OBLIQUE incidence on the output of
crack_oblique.mac: 16 fixed-gun runs at 45 deg off normal (IWaveTheta = 135)
into crack1, one ROOT file per run (output/bbr_oblique_r00.root ... r15).

Physics the checks rest on (independent of the code under test):
  * The parallel-plate gap is translation-invariant along its long axis (y),
    so a transmitted photon keeps the SIGN of its incoming k_y. The HFSS
    far-field tables carry this (the mean outgoing y-component tracks the
    incoming one). Normal incidence cannot see it because k_y = 0.
  * The gap is a polarization filter: E across the gap (z) is the TEM
    polarization and transmits; at a y-tilt, E in the plane of incidence has
    no z-component and is cut off. At a z-tilt the roles swap.
  * Mirror symmetry: +y/-y tilts (and +z/-z) give mirror-image exit patterns.
  * Reflection is specular about the crack normal (x).
Quantitative predictions (T, mean exit direction, exit marginals) come from
bbrsim.hfss, a Python mirror of BBRHFSSData reading the same CSVs.

Per run: exactly one crack1 entry per event; entry direction equals the
configured gun direction; T_obs within 3 sigma of the prediction (fixed
polarization: the transmitting or cut-off value from the filter physics);
every reflected photon exactly specular. Random-polarization runs add: mean
exit direction within 4 SE of the mixture prediction per component;
sign(mean k_y,out) = sign(k_y,in); chi^2 of the k_y and k_z exit marginals
against the predicted table (p > 1e-3); two-sample KS mirror tests between
partner runs (p > 0.01).

Usage:
    conda run -n bbrsim python scripts/check_crack_oblique.py [output_dir | glob]
    (default: build/output)
"""
import glob
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
from bbrsim import hfss, select
from bbrsim.io import load_many

DIR_STEM = "InfParallelPlate_crack1Rohan_500GHz"
S = 0.70710678

# Must match the run order in crack_oblique.mac: (label, gun dir, pol, events).
# pol None = random polarization; otherwise the requested vector, which the
# PGA projects perpendicular to the direction and normalises.
RUNS = [
    ("+y   random",         (S,  S, 0.), None,       8000),
    ("-y   random",         (S, -S, 0.), None,       8000),
    ("+z   random",         (S, 0.,  S), None,       8000),
    ("-z   random",         (S, 0., -S), None,       8000),
    ("+y+z random",         (S, .5, .5), None,       8000),
    ("+y-z random",         (S, .5, -.5), None,      8000),
    ("-y+z random",         (S, -.5, .5), None,      8000),
    ("-y-z random",         (S, -.5, -.5), None,     8000),
    ("+y   E||z  (TEM)",    (S,  S, 0.), (0., 0., 1.), 3000),
    ("+y   E in-plane",     (S,  S, 0.), (0., 1., 0.), 3000),
    ("-y   E||z  (TEM)",    (S, -S, 0.), (0., 0., 1.), 3000),
    ("-y   E in-plane",     (S, -S, 0.), (0., 1., 0.), 3000),
    ("+z   E||y  (cutoff)", (S, 0.,  S), (0., 1., 0.), 3000),
    ("+z   E in-plane",     (S, 0.,  S), (0., 0., 1.), 3000),
    ("-z   E||y  (cutoff)", (S, 0., -S), (0., 1., 0.), 3000),
    ("-z   E in-plane",     (S, 0., -S), (0., 0., 1.), 3000),
]
# (run a, run b, mirrored axis): b is a's mirror image in that axis.
MIRROR_PAIRS = [(0, 1, "y"), (2, 3, "z"), (4, 6, "y"), (4, 5, "z"), (7, 5, "y"), (7, 6, "z")]

Z_MEAN_MAX = 4.0     # SE units, mean exit direction per component
SIGMA_T_MAX = 3.0    # binomial sigma, transmittance
P_CHI2_MIN = 1e-3    # exit-marginal chi^2
P_KS_MIN = 0.01      # mirror KS
TOL_DIR = 1e-6       # entry direction vs configured
TOL_REFL = 1e-9      # specular reflection

arg = sys.argv[1] if len(sys.argv) > 1 else "build/output"
files = sorted(glob.glob(os.path.join(arg, "bbr_oblique_r*.root")) if os.path.isdir(arg)
               else glob.glob(arg))
if not files:
    print(f"no ROOT files matched {arg!r}")
    print("RESULT: FAIL")
    sys.exit(1)

df = load_many(files)
entries = select.crack_crossings(df)
entries = entries[entries["vol_post"].str.contains("crack1", na=False)]
datasets = hfss.load_dataset(DIR_STEM)

results = []          # (name, passed, detail)
per_run = {}          # run index -> dict for mirror tests


def check(name, passed, detail=""):
    results.append((name, bool(passed), detail))


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def chi2_marginal(values, pred_values, pred_weights, n_bins=12, decimals=6):
    """Chi^2 of observed ``values`` against the predicted weighted table.

    The far-field grid is coarse (37 Theta values), so a direction component
    such as k_z = cos(Theta) is a DISCRETE distribution with large atoms. Bin
    edges are therefore placed midway BETWEEN predicted atoms, never on one,
    and both samples are rounded to ``decimals`` so a last-bit difference
    between the C++ and numpy cosines cannot move an atom across an edge.
    """
    pv = np.round(pred_values, decimals)
    ov = np.round(values, decimals)
    atoms, inv = np.unique(pv, return_inverse=True)
    mass = np.bincount(inv, weights=pred_weights)
    cum = np.cumsum(mass)
    cuts = np.searchsorted(cum, np.linspace(0., 1., n_bins + 1)[1:-1])
    cuts = np.unique(cuts[(cuts >= 0) & (cuts < len(atoms) - 1)])
    edges = np.concatenate(([-np.inf], 0.5 * (atoms[cuts] + atoms[cuts + 1]), [np.inf]))
    if len(edges) < 4:
        return float("nan"), 0
    expected = np.array([mass[(atoms > lo) & (atoms <= hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])])
    expected *= len(ov)
    observed = np.array([((ov > lo) & (ov <= hi)).sum() for lo, hi in zip(edges[:-1], edges[1:])])
    keep = expected > 0
    chi2 = float(((observed[keep] - expected[keep]) ** 2 / expected[keep]).sum())
    ndf = int(keep.sum()) - 1
    return stats.chi2.sf(chi2, ndf) if ndf > 0 else float("nan"), ndf


print(f"files: {len(files)}   crossings: {len(df)}   crack1 entries: {len(entries)}")
print(f"{'run':>3} {'config':<20} {'N':>5} {'T_obs':>7} {'T_pred':>7} {'z_T':>5} "
      f"{'<py>obs':>8} {'<py>pred':>8} {'<pz>obs':>8}")

for r, (label, gdir, pol, n_events) in enumerate(RUNS):
    rows = entries[entries["run_id"] == r]
    N = len(rows)
    tag = f"run{r:02d}"
    if N == 0:
        check(f"{tag} present", False, "no crack1 entries")
        continue

    # --- geometry self-checks ---------------------------------------------
    one_each = (rows.groupby("event_id").size() == 1).all() and N == n_events
    check(f"{tag} one entry per event", one_each, f"{N} entries for {n_events} events")

    k_in = rows[["px_pre", "py_pre", "pz_pre"]].to_numpy(float)
    k_cfg = unit(gdir)
    spread = np.abs(k_in - k_cfg).max()
    check(f"{tag} entry direction", spread < TOL_DIR, f"max |k - k_cfg| = {spread:.1e}")
    k0 = k_in[0]

    # --- predictions from the Python mirror --------------------------------
    inc = hfss.fold_incidence(k0)
    key = hfss.nearest_key(datasets, inc.phi_deg, inc.theta_deg)
    ds = datasets[key]
    dirs_pred = hfss.outgoing_directions(ds, inc)
    T_tem, T_cut = max(ds.T0, ds.T1), min(ds.T0, ds.T1)

    if pol is None:
        T_pred, w_pred = hfss.random_polarization_mixture(ds)
        pol_note = "random"
    else:
        p = np.asarray(pol, float)
        p_perp = p - np.dot(p, k0) * k0          # what the PGA does
        p_perp = unit(p_perp)
        # Filter physics: only E with a component across the gap (z) transmits.
        T_pred = T_tem if abs(p_perp[2]) > 0.5 else T_cut
        Et, Ep = hfss.polarization_components(p_perp, inc)
        T_basis = hfss.transmittance(ds, Et, Ep)   # what the code's basis gives
        w_pred = hfss.direction_weights(ds, Et, Ep)
        pol_note = f"E={np.round(p_perp, 3)} basis-T={T_basis:.4f}"
        check(f"{tag} basis agrees with filter physics", abs(T_basis - T_pred) < 1e-3,
              f"basis {T_basis:.4f} vs physics {T_pred:.4f}")

    # --- transmittance ------------------------------------------------------
    tx = rows[rows["status"] == "BBRDiffractionTransmit"]
    rf = rows[rows["status"] == "BBRDiffractionReflect"]
    k_tx = len(tx)
    T_obs = k_tx / N
    if T_pred < 1e-6:
        z_T = float("nan")
        check(f"{tag} transmittance", k_tx <= 1, f"{k_tx} transmitted, expected ~0 (cut-off)")
    else:
        sigma = np.sqrt(T_pred * (1 - T_pred) / N)
        z_T = (T_obs - T_pred) / sigma
        check(f"{tag} transmittance", abs(z_T) < SIGMA_T_MAX,
              f"T_obs {T_obs:.4f} vs {T_pred:.4f} ({z_T:+.2f} sigma) [{pol_note}]")
    check(f"{tag} all entries decided", k_tx + len(rf) == N,
          f"{k_tx} transmitted + {len(rf)} reflected, {N} entries")

    # --- reflection is specular about x --------------------------------------
    if len(rf):
        d_post = rf[["px_post", "py_post", "pz_post"]].to_numpy(float)
        d_exp = np.array([hfss.reflected_direction(k) for k in rf[["px_pre", "py_pre", "pz_pre"]].to_numpy(float)])
        dev = np.abs(d_post - d_exp).max()
        check(f"{tag} reflection specular", dev < TOL_REFL, f"max deviation {dev:.1e}")

    # --- exit direction (transmitted photons) --------------------------------
    d_out = tx[["px_post", "py_post", "pz_post"]].to_numpy(float) if k_tx else np.zeros((0, 3))
    mean_pred = w_pred @ dirs_pred
    mean_obs = d_out.mean(axis=0) if k_tx else np.full(3, np.nan)
    per_run[r] = {"py": d_out[:, 1] if k_tx else np.array([]),
                  "pz": d_out[:, 2] if k_tx else np.array([])}

    if pol is None and k_tx > 100:
        se = d_out.std(axis=0, ddof=1) / np.sqrt(k_tx)
        z = (mean_obs - mean_pred) / np.where(se > 0, se, np.inf)
        check(f"{tag} mean exit direction", np.all(np.abs(z) < Z_MEAN_MAX),
              f"obs {np.round(mean_obs, 4)} pred {np.round(mean_pred, 4)} z {np.round(z, 2)}")
        # k_y sign conservation along the invariant (long) axis
        if abs(k0[1]) > 0.1:
            check(f"{tag} k_y sign conserved", np.sign(mean_obs[1]) == np.sign(k0[1]),
                  f"k_y,in {k0[1]:+.3f}  <k_y,out> {mean_obs[1]:+.4f}")
        else:
            check(f"{tag} <k_y,out> ~ 0", abs(mean_obs[1]) < Z_MEAN_MAX * se[1],
                  f"<k_y,out> {mean_obs[1]:+.4f} (SE {se[1]:.4f})")
        if abs(k0[2]) > 0.1:
            check(f"{tag} <k_z,out> ~ 0 (symmetric gap)", abs(mean_obs[2]) < Z_MEAN_MAX * se[2],
                  f"<k_z,out> {mean_obs[2]:+.4f} (SE {se[2]:.4f})")
        for comp, ci in (("k_y", 1), ("k_z", 2)):
            p, ndf = chi2_marginal(d_out[:, ci], dirs_pred[:, ci], w_pred)
            check(f"{tag} {comp} exit marginal chi2", (not np.isnan(p)) and p > P_CHI2_MIN,
                  f"p = {p:.3g} (ndf {ndf})")

    print(f"{r:>3} {label:<20} {N:>5} {T_obs:>7.4f} {T_pred:>7.4f} {z_T:>+5.1f} "
          f"{mean_obs[1]:>+8.4f} {mean_pred[1]:>+8.4f} {mean_obs[2]:>+8.4f}")

# --- mirror symmetry between partner runs ------------------------------------
for a, b, axis in MIRROR_PAIRS:
    if a not in per_run or b not in per_run:
        continue
    A, B = per_run[a], per_run[b]
    if len(A["py"]) < 100 or len(B["py"]) < 100:
        continue
    flip, same = ("py", "pz") if axis == "y" else ("pz", "py")
    p_flip = stats.ks_2samp(A[flip], -B[flip]).pvalue
    p_same = stats.ks_2samp(A[same], B[same]).pvalue
    check(f"mirror run{a:02d}/run{b:02d} in {axis}: {flip} flips", p_flip > P_KS_MIN, f"KS p = {p_flip:.3g}")
    check(f"mirror run{a:02d}/run{b:02d} in {axis}: {same} same", p_same > P_KS_MIN, f"KS p = {p_same:.3g}")

# --- report --------------------------------------------------------------------
print()
n_fail = 0
for name, ok, detail in results:
    if not ok:
        n_fail += 1
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<48} {detail}")
print(f"\nchecks: {len(results)}   failed: {n_fail}")
print(f"RESULT: {'PASS' if n_fail == 0 and results else 'FAIL'}")
sys.exit(0 if n_fail == 0 and results else 1)
