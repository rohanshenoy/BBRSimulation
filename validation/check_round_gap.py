"""
check_round_gap.py
Validate the straight round gap on the output of Validation_RoundGap.mac.

Setup: a vacuum_wg G4DisplacedSolid tube (radius 51 um, axis along x) in a
0.4 mm Cu plate at z = -80 mm, fed by the mock HFSS dataset
RoundGap_r50um_2000GHz (validation/Scripts/make_mock_round_gap.py). Nine
fixed-gun runs at 2000 GHz write output/bbr_round_r00.root ... r08: normal
incidence with E across g, E along l and random polarization; the four 45 deg
tilts toward +-l and +-g; the diagonal (IWavePhi 45, IWaveTheta 135) with random
polarization, where the cross term averages out; and the same diagonal with a
fixed polarization of equal theta-hat and phi-hat components, where the mock's
Re rho = 0.288 makes the cross term 2 Et Ep sqrt(T0 T1) Re rho large. There
T = 0.511, and 0.289 with the cross term's sign flipped, 28 sigma apart at 4000
events, so r08 pins that sign end to end.

Per run:
  * one RoundGap_r50um entry per event, along the configured direction;
  * T_obs within 3 sigma of the bbrsim.hfss prediction (fixed or random polarization);
  * every exit on the exit face (x = 0.4 mm) and within the HFSS radius 50 um of
    the axis: positions come from the table, not from the 51 um solid;
  * fixed polarization only: the share of exits with r < R/2 is within 4 sigma
    of the table's prediction. The Ephi=0 profile (1 - (r/R)^2)^2 and the Ephi=1
    profile (r/R)^4 differ, so this sees the polarization mapping;
  * the mean exit direction is within 4 SE of the prediction in each component.
With --log <run.log>: the startup check (BBRCrackLibrary::ValidatePlacedCracks)
printed exactly one "[BBR] crack RoundGap_r50um:" line reporting its sidecar(s)
fit, so the placed solid was checked against the sidecar before any photon.

49 checks with --log, 48 without: 5 for each of the 6 random-polarization runs
and 6 for each of the 3 fixed-polarization runs.

Usage:
    conda run -n bbrsim python validation/check_round_gap.py <output_dir> --data-dir <mock root> [--log run.log]
"""
import argparse
import os
import re
import sys

import numpy as np

from bbrsim import hfss
from bbrsim.io import load

ID, STEM = "RoundGap_r50um", "RoundGap_r50um_2000GHz"
R_MM, Z0_MM, EXIT_X_MM = 0.050, -80.0, 0.4
S = 0.70710678
RUNS = [   # (label, gun direction, polarization or None = random, events): the macro's order
    ("normal E||g (Ephi=0)", (1., 0., 0.), (0., 0., 1.), 4000),
    ("normal E||l (Ephi=1)", (1., 0., 0.), (0., 1., 0.), 4000),
    ("normal random", (1., 0., 0.), None, 4000),
    ("+l tilt random", (S, S, 0.), None, 4000),
    ("-l tilt random", (S, -S, 0.), None, 4000),
    ("+g tilt random", (S, 0., S), None, 4000),
    ("-g tilt random", (S, 0., -S), None, 4000),
    ("diagonal random", (S, .5, .5), None, 4000),
    ("diagonal fixed (Et = Ep)", (S, .5, .5), (.5, .14644661, -.85355339), 4000),
]
SIGMA_T_MAX, Z_MAX = 3.0, 4.0

ap = argparse.ArgumentParser()
ap.add_argument("output", nargs="?", default="output")
ap.add_argument("--data-dir", required=True, help="data root holding waveguides/RoundGap_r50um_2000GHz_*")
ap.add_argument("--log", help="the run log; check the startup crack line of RoundGap_r50um")
args = ap.parse_args()
base = os.path.join(args.data_dir, "waveguides")
datasets = hfss.load_dataset(STEM, base)
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def inner_share(ds, Et, Ep):
    """Predicted share of exits with r < R/2 for (E_theta, E_phi) at one key."""
    r_mm = np.hypot(ds.exit_y_m, ds.exit_z_m) * 1e3
    return float(hfss.exit_position_weights(ds, Et, Ep)[r_mm < R_MM / 2].sum())


if args.log:
    try:
        with open(args.log, errors="replace") as fh:
            lines = [ln.rstrip() for ln in fh if f"[BBR] crack {ID}:" in ln]
    except OSError as err:
        lines = [f"unreadable: {err}"]
    fit = [re.search(r"; (\d+) sidecar\(s\) fit$", ln) for ln in lines]
    check("startup crack check in the log", len(lines) == 1 and fit[0] is not None and int(fit[0].group(1)) >= 1,
          lines[0] if len(lines) == 1 else f"{len(lines)} '[BBR] crack {ID}:' lines in {args.log}")


for i, (label, k, pol, n) in enumerate(RUNS):
    tag = f"r{i:02d} {label}"
    path = os.path.join(args.output, f"bbr_round_r{i:02d}.root")
    if not os.path.exists(path):
        check(f"{tag}: file", False, path)
        continue
    cr, _ = load(path)
    ent = cr[(cr.mat_post == "vacuum_wg") & (cr.vol_post == ID)]
    check(f"{tag}: one entry per event", len(ent) == n and ent.event_id.is_unique, f"{len(ent)} entries, {n} events")
    kin = np.stack([ent.px_pre, ent.py_pre, ent.pz_pre], axis=1)
    check(f"{tag}: entry direction", np.allclose(kin, unit(k), atol=1e-6))
    inc = hfss.fold_incidence(unit(k))
    ds = datasets[hfss.nearest_key(datasets, inc.phi_deg, inc.theta_deg)]
    share = None
    if pol is None:
        T_pred, W = hfss.random_polarization_mixture(ds)
    else:
        p = np.asarray(pol, float)
        p = p - np.dot(p, unit(k)) * unit(k)
        Et, Ep = hfss.polarization_components(unit(p), inc)
        T_pred, W = hfss.transmittance(ds, Et, Ep), hfss.direction_weights(ds, Et, Ep)
        share = inner_share(ds, Et, Ep)
    tr = ent[ent.status == "BBRDiffractionTransmit"]
    T_obs = len(tr) / max(len(ent), 1)
    sig = np.sqrt(max(T_pred * (1 - T_pred), 1e-12) / max(len(ent), 1))
    check(f"{tag}: T", abs(T_obs - T_pred) <= SIGMA_T_MAX * sig, f"T_obs {T_obs:.4f} vs {T_pred:.4f} (sigma {sig:.4f})")
    ex = cr[(cr.mat_pre == "vacuum_wg") & (cr.vol_pre == ID) & cr.event_id.isin(tr.event_id)]
    r = np.hypot(ex.y_mm.to_numpy(), ex.z_mm.to_numpy() - Z0_MM)
    check(f"{tag}: exits on the face, inside R", len(ex) == len(tr) and bool(np.all(r <= R_MM * (1 + 1e-6)))
          and bool(np.allclose(ex.x_mm, EXIT_X_MM, atol=1e-6)), f"{len(ex)} exits, max r {r.max() if len(r) else 0:.6f} mm")
    if share is not None and len(r):
        obs = float(np.mean(r < R_MM / 2))
        s2 = np.sqrt(max(share * (1 - share), 1e-12) / len(r))
        check(f"{tag}: radial share r < R/2", abs(obs - share) <= Z_MAX * s2, f"{obs:.4f} vs {share:.4f}")
    D = hfss.outgoing_directions(ds, inc)
    mean_pred = W @ D
    sd = np.sqrt(np.maximum(W @ (D - mean_pred) ** 2, 1e-30))
    kout = np.stack([tr.px_post, tr.py_post, tr.pz_post], axis=1)
    if len(kout):
        z = np.abs(kout.mean(axis=0) - mean_pred) / (sd / np.sqrt(len(kout)))
        check(f"{tag}: mean exit direction", bool(np.all(z <= Z_MAX)),
              f"obs {np.round(kout.mean(axis=0), 4)} pred {np.round(mean_pred, 4)} z {np.round(z, 2)}")

n_fail = sum(1 for _, ok, _ in results if not ok)
for name, ok, detail in results:
    print(f"  {'PASS' if ok else 'FAIL'} {name}  {detail}")
print(f"{len(results)} checks, {n_fail} failing")
print("RESULT: " + ("PASS" if results and n_fail == 0 else "FAIL"))
sys.exit(0 if results and n_fail == 0 else 1)
