#!/usr/bin/env python3
"""
make_tree_check.py: write a test-world run that exercises every HFSS dataset of a
data tree the test world can place, for validation/check_tree_run.py.

    conda run -n bbrsim python validation/Scripts/make_tree_check.py DATA_ROOT --out RUN_DIR
        [--events 4000] [--probe-events 1000] [--planck-events 100000]
        [--planck-temp-K 20] [--seeds 2718 28]

DATA_ROOT holds waveguides/. Only directory names and sidecars are read, never
the CSVs, so check the tree first with validation/check_dataset_sidecars.py.
Writes RUN_DIR/tree_check.mac and RUN_DIR/tree_check.json (the run manifest);
then, from RUN_DIR:

    G4FORCENUMBEROFTHREADS=8 bbrsimTestWorld tree_check.mac > run.log 2>&1
    conda run -n bbrsim python validation/check_tree_run.py RUN_DIR --log RUN_DIR/run.log

The run plan, per dataset the test world places (InfParallelPlate_crack1Rohan,
InfParallelPlate_crack2, and RoundGap_r50um through /bbr/testworld/roundGap)
and per grid frequency, gun photons from x = -20 mm aimed at the entrance-face
centre:
  normal_random, normal_E_gap, normal_E_long: normal incidence with random
      polarization, E across the gap, E along the long side;
  oblique_random: every declared IWaveTheta strictly between 90 and 180 deg
      and every declared IWavePhi, random polarization;
  oblique_diagonal_fixed: IWavePhi 45 on each such row, E_theta = E_phi.
On grids of more than one point the lowest and highest grid frequency run
0.1 % inside the grid, and probes go to 0.99 and 1.01 times each log midpoint
and to half the lowest and twice the highest grid frequency. Each dataset gets
one Planck run from a thin box just in front of its crack. Datasets the test
world cannot place are listed as not exercised. Keep every event count at or
above the square of the Geant4 thread count (Geant4 warns, Run10035, otherwise).
Exit code 0 on success, 2 when the tree cannot be checked.
"""
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys

import numpy as np

from bbrsim import hfss, sidecar

SCHEMA = "bbrsim-tree-check/1"
GUN_X_MM = -20.0
EMITTER_X_MM = -0.03
PLANCK_BAND_EV = (4.14e-5, 8.27e-2)        # BBRPrimarySource's fixed emitter band
EDGE_NUDGE = 1e-3                          # grid-edge runs sit 0.1 % inside the grid
PROBE_FACTOR = 0.01                        # probes at (1 -/+ 0.01) x each log midpoint
NORMAL_KEY = (0.0, 180.0)                  # (IWavePhi, IWaveTheta)
E_GAP, E_LONG = (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)
# The test world's cracks (examples/testworld/src/TestWorldDetectorConstruction.cc):
# unrotated, propagation along world +x, long side along +y, gap along +z.
AIM = {
    "InfParallelPlate_crack1Rohan": {"face_mm": (0.0, 0.0, 0.0), "emitter_mm": (0.02, 8.0, 0.04)},
    "InfParallelPlate_crack2": {"face_mm": (0.0, 0.0, 3.0), "emitter_mm": (0.02, 8.0, 0.08)},
    "RoundGap_r50um": {"face_mm": (0.0, 0.0, -80.0), "emitter_mm": (0.02, 0.07, 0.07),
                       "switch": "/bbr/testworld/roundGap true"},
}
ALWAYS_PLACED = ("InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2")


class Refusal(Exception):
    """The tree cannot be checked; the message names the problem and the fix."""


def direction(theta_deg, phi_deg):
    """World unit vector of a photon arriving at (IWavePhi, IWaveTheta) on an
    unrotated test-world crack: 180 - theta off +x, azimuth phi from the gap axis
    (+z) towards the long axis (+y), so that it folds back onto (phi, theta)."""
    a, p = math.radians(180.0 - theta_deg), math.radians(phi_deg)
    k = (math.cos(a), math.sin(a) * math.sin(p), math.sin(a) * math.cos(p))
    return tuple(0.0 if abs(c) < 1e-12 else c for c in k)


def start_mm(face_mm, k):
    """Gun position on the line through the entrance-face centre, at x = GUN_X_MM."""
    t = (face_mm[0] - GUN_X_MM) / k[0]
    return tuple(f - t * c for f, c in zip(face_mm, k))


def diagonal_polarization(k):
    """(e_theta + e_phi) / sqrt(2) of the wrapper's incident basis at k: E_theta = E_phi."""
    e_t, e_p = hfss.incoming_basis(hfss.fold_incidence(k))
    return tuple(float(c) for c in (np.asarray(e_t) + np.asarray(e_p)) / math.sqrt(2.0))


def energy_eV(f_ghz):
    """Gun energy for frequency f_ghz, rounded to the 11 significant figures the macro carries."""
    return float(f"{f_ghz * 1e9 * hfss.H_EV_S:.10e}")


def read_tree(root):
    """(waveguides dir, {id: manifest entry}) for the tree under root; raises Refusal."""
    wg = os.path.join(root, "waveguides")
    if not os.path.isdir(wg):
        raise Refusal(f"{wg} is not a directory: DATA_ROOT must hold waveguides/")
    try:
        grids = hfss.discover_datasets(wg)
    except ValueError as err:
        raise Refusal(f"{wg}: {err}") from None
    missing = [i for i in ALWAYS_PLACED if i not in grids]
    if missing:
        raise Refusal(f"{wg} has no dataset for {', '.join(missing)}: the test world always places "
                      "both slab cracks and checks every placed crack at /run/initialize")
    datasets = {}
    for did, grid in grids.items():
        entry = {"aimed": did in AIM, "grid": [[f, stem] for f, stem in grid],
                 "sidecar_sha256": {stem: sidecar_sha256(wg, stem) for _, stem in grid}}
        if did in AIM:
            try:
                sc = [sidecar.load(wg, stem) for _, stem in grid][0]
                # the frequency-independent blocks agree across a grid (check_dataset_sidecars.py)
                entry.update({
                    "cutoff_ghz": float(sc["modes"]["cutoff_ghz"]),
                    "polarization_filter_limit_ghz": float(sc["modes"]["polarization_filter_limit_ghz"]),
                    "cross_section": str(sc["exit_field"]["cross_section"]["shape"]),
                    "incident_phi_deg": [float(x) for x in sc["excitation"]["incident_phi_deg"]],
                    "incident_theta_deg": [float(x) for x in sc["excitation"]["incident_theta_deg"]]})
            except ValueError as err:          # BBR024 from sidecar.load, or a non-number
                raise Refusal(f"{err}; check the tree with validation/check_dataset_sidecars.py {wg}") from None
            except (KeyError, TypeError) as err:
                raise Refusal(f"{sidecar.path_for(wg, grid[0][1])}: field {err} missing or malformed; "
                              f"check the tree with validation/check_dataset_sidecars.py {wg}") from None
        datasets[did] = entry
    return wg, datasets


def sidecar_sha256(wg, stem):
    """sha256 of a dataset's sidecar, which records its CSVs' checksums; None if unreadable."""
    try:
        with open(os.path.join(wg, f"{stem}.dataset.json"), "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return None


def plan_runs(datasets, events, probe_events, planck_events, planck_temp_K):
    """(runs, skipped): the manifest's runs in macro order, and the runs not generated."""
    runs, skipped = [], []

    def gun(kind, did, f_target, k, pol, n, key):
        inc = hfss.fold_incidence(k)
        if abs(inc.phi_deg - key[0]) > 1e-6 or abs(inc.theta_deg - key[1]) > 1e-6:
            raise RuntimeError(f"direction {k} folds to ({inc.phi_deg}, {inc.theta_deg}), not {key}")
        e = energy_eV(f_target)
        grid = [(f, s) for f, s in datasets[did]["grid"]]
        i = len(runs)
        runs.append({"run": i, "file": f"output/tree_r{i:03d}.root", "kind": kind, "dataset": did,
                     "target_ghz": f_target, "energy_eV": e,
                     "expected_grid_ghz": hfss.select_frequency(grid, hfss.photon_frequency_GHz(e))[0],
                     "events": n, "position_mm": list(start_mm(AIM[did]["face_mm"], k)),
                     "direction": list(k), "polarization": None if pol is None else list(pol),
                     "incidence_key": list(key)})

    aimed = sorted(d for d in datasets if datasets[d]["aimed"])
    normal = direction(180.0, 0.0)
    for did in aimed:
        ds = datasets[did]
        freqs = [f for f, _ in ds["grid"]]
        rows = [t for t in ds["incident_theta_deg"] if 90.0 < t < 180.0]
        phis = [p for p in ds["incident_phi_deg"] if 0.0 <= p <= 90.0]
        if not rows:
            skipped.append({"dataset": did, "kinds": "oblique",
                            "reason": "no declared IWaveTheta strictly between 90 and 180 deg"})
        for t in rows:
            if 45.0 not in phis:
                skipped.append({"dataset": did, "kinds": "oblique_diagonal_fixed",
                                "reason": f"IWavePhi 45 not declared (IWaveTheta {t:g})"})
        multi = len(freqs) > 1
        for i, f in enumerate(freqs):
            ft = f
            if multi and i == 0:
                ft = f * (1 + EDGE_NUDGE)
            elif multi and i == len(freqs) - 1:
                ft = f * (1 - EDGE_NUDGE)
            gun("normal_random", did, ft, normal, None, events, NORMAL_KEY)
            gun("normal_E_gap", did, ft, normal, E_GAP, events, NORMAL_KEY)
            gun("normal_E_long", did, ft, normal, E_LONG, events, NORMAL_KEY)
            for t in rows:
                for p in phis:
                    gun("oblique_random", did, ft, direction(t, p), None, events, (p, t))
                if 45.0 in phis:
                    k = direction(t, 45.0)
                    gun("oblique_diagonal_fixed", did, ft, k, diagonal_polarization(k), events, (45.0, t))
        if multi:
            for lo, hi in zip(freqs[:-1], freqs[1:]):
                m = math.sqrt(lo * hi)
                for fac in (1 - PROBE_FACTOR, 1 + PROBE_FACTOR):
                    gun("probe_midpoint", did, m * fac, normal, None, probe_events, NORMAL_KEY)
            gun("probe_clamp_low", did, freqs[0] / 2.0, normal, None, probe_events, NORMAL_KEY)
            gun("probe_clamp_high", did, freqs[-1] * 2.0, normal, None, probe_events, NORMAL_KEY)
    for did in aimed:
        face, i = AIM[did]["face_mm"], len(runs)
        runs.append({"run": i, "file": f"output/tree_r{i:03d}.root", "kind": "planck", "dataset": did,
                     "events": planck_events, "temperature_K": planck_temp_K,
                     "emitter_center_mm": [EMITTER_X_MM, face[1], face[2]],
                     "emitter_size_mm": list(AIM[did]["emitter_mm"])})
    return runs, skipped


def fmt(x):
    """Shortest round-trip decimal; the /bbr/gun and /bbr/thermal commands parse it in full."""
    return repr(float(x))


def macro_text(root, datasets, runs, seeds):
    lines = ["# tree_check.mac, written by validation/Scripts/make_tree_check.py; do not edit.",
             f"# Data root {root}; {len(runs)} runs, listed in tree_check.json.",
             "/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0",
             f'/bbr/dataDir "{root}"']
    lines += [AIM[d]["switch"] for d in sorted(datasets) if datasets[d]["aimed"] and "switch" in AIM[d]]
    lines += [f"/random/setSeeds {seeds[0]} {seeds[1]}", "/run/initialize"]
    for r in runs:
        what = f"{r['target_ghz']:.6g} GHz" if "target_ghz" in r else f"{r['temperature_K']:g} K"
        lines.append(f"# run {r['run']}: {r['dataset']} {r['kind']} {what}")
        if r["kind"] == "planck":
            c, s = r["emitter_center_mm"], r["emitter_size_mm"]
            lines += ["/bbr/gun/mode false", f"/bbr/thermal/setT {fmt(r['temperature_K'])} K",
                      f"/bbr/thermal/emitterCenter {fmt(c[0])} {fmt(c[1])} {fmt(c[2])} mm",
                      f"/bbr/thermal/emitterSize {fmt(s[0])} {fmt(s[1])} {fmt(s[2])} mm"]
        else:
            p, k = r["position_mm"], r["direction"]
            pol = r["polarization"] or (0.0, 0.0, 0.0)
            lines += ["/bbr/gun/mode true", f"/bbr/gun/energy_eV {r['energy_eV']:.10e}",
                      f"/bbr/gun/posX {fmt(p[0])}", f"/bbr/gun/posY {fmt(p[1])}", f"/bbr/gun/posZ {fmt(p[2])}",
                      f"/bbr/gun/dirX {fmt(k[0])}", f"/bbr/gun/dirY {fmt(k[1])}", f"/bbr/gun/dirZ {fmt(k[2])}",
                      f"/bbr/gun/pol {fmt(pol[0])} {fmt(pol[1])} {fmt(pol[2])}"]
        lines += [f"/analysis/setFileName {r['file'][:-len('.root')]}", f"/run/beamOn {r['events']}"]
    return "\n".join(lines) + "\n"


def git_describe():
    try:
        r = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "describe",
                            "--always", "--dirty"], capture_output=True, text=True, timeout=10)
        return r.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def write_files(out, files):
    """Write every (name, text) under a temporary name first, then rename them all."""
    os.makedirs(out, exist_ok=True)
    staged = []
    for name, text in files:
        tmp = os.path.join(out, f".{name}.tmp")
        with open(tmp, "w") as fh:
            fh.write(text)
        staged.append((tmp, os.path.join(out, name)))
    for tmp, final in staged:
        os.replace(tmp, final)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_root", help="data root holding waveguides/")
    ap.add_argument("--out", required=True, help="run directory for tree_check.mac and tree_check.json")
    ap.add_argument("--events", type=int, default=4000, help="photons per gun run")
    ap.add_argument("--probe-events", type=int, default=1000, help="photons per frequency probe")
    ap.add_argument("--planck-events", type=int, default=100000, help="photons per Planck run")
    ap.add_argument("--planck-temp-K", type=float, default=20.0, help="Planck emitter temperature")
    ap.add_argument("--seeds", type=int, nargs=2, default=(2718, 28), help="/random/setSeeds")
    a = ap.parse_args(argv)
    root, out = os.path.realpath(a.data_root), os.path.realpath(a.out)
    try:
        if out == root or out.startswith(root + os.sep):
            raise Refusal(f"refusing to write into the data tree: {out} lies inside {root}")
        if min(a.events, a.probe_events, a.planck_events) <= 0 or not a.planck_temp_K > 0:
            raise Refusal("event counts and the Planck temperature must be positive")
        wg, datasets = read_tree(root)
        runs, skipped = plan_runs(datasets, a.events, a.probe_events, a.planck_events, a.planck_temp_K)
    except Refusal as err:
        print(f"make_tree_check: {err}", file=sys.stderr)
        return 2
    not_exercised = sorted(d for d in datasets if not datasets[d]["aimed"])
    manifest = {"schema": SCHEMA, "data_root": root, "waveguides": wg, "seeds": list(a.seeds),
                "h_eV_s": hfss.H_EV_S,
                "planck": {"temperature_K": a.planck_temp_K, "band_eV": list(PLANCK_BAND_EV)},
                "generator_version": git_describe(), "datasets": datasets,
                "not_exercised": not_exercised, "skipped": skipped, "runs": runs}
    write_files(out, [("tree_check.mac", macro_text(root, datasets, runs, a.seeds)),
                      ("tree_check.json", json.dumps(manifest, indent=1) + "\n")])
    n_aimed = sum(d["aimed"] for d in datasets.values())
    print(f"wrote {len(runs)} runs for {n_aimed} dataset(s) to {out}"
          + (f"; not exercised: {', '.join(not_exercised)}" if not_exercised else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
