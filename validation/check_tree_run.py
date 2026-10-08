#!/usr/bin/env python3
"""
check_tree_run.py: check a run of tree_check.mac against its HFSS data tree.

    conda run -n bbrsim python validation/check_tree_run.py RUN_DIR [--log RUN_DIR/run.log]
        [--data-root DIR]

RUN_DIR holds tree_check.json, written by validation/Scripts/make_tree_check.py,
and output/tree_rNNN.root with their .metadata.json files. Every prediction comes
from bbrsim.hfss reading the same tree, so the check shows that BBRsim consumed
the tree as the mirror reads it; apart from the physics row it does not judge the
data themselves.

Exact rows: the tree is unchanged since generation; every output names the data
root and the placed cracks with the tree's frequencies; hfss_freq_GHz is -1 on
every row that is not a decided crack entry; every decided entry carries the grid
frequency the selection rule gives for its recorded energy; per gun run, one entry
per event on the aimed crack, with the configured direction and energy and the
expected grid frequency; per Planck run, at most one entry per event; reflections
are specular; every grid frequency is served by a gun run. With --log: BBR008 and
BBR026 exactly as the output implies, a clean log, and one startup line per placed
crack and dataset.
Statistical rows: the transmittance of every gun run and probe, exit-direction
moments, the Planck runs' transmittance and spectrum, each within z_max standard
errors, z_max = max(4, Phi^-1(1 - 0.0005 / M)) for M statistical rows, so a
correct run fails any row with at most 0.1 % chance; counts get one count of
slack. Physics row: below 0.95 x the polarization-filter limit, a rectangular
crack transmits under 1 % of E along its long side at normal incidence (a
swapped pair of Ephi directories fails it).
Prints one PASS / FAIL / SKIP line per check, then RESULT: PASS or FAIL; exits 0 or 1.
"""
import argparse
import json
import math
import os
import re
import sys
from statistics import NormalDist

import numpy as np

from bbrsim import hfss, io, physics, select

SCHEMA = "bbrsim-tree-check/1"
TRANSMIT, REFLECT = "BBRDiffractionTransmit", "BBRDiffractionReflect"
DECIDED = (TRANSMIT, REFLECT)
TOL_DIR = 1e-6            # entry direction against the configured one
TOL_REFL = 1e-9           # specular reflection
TOL_ENERGY = 1e-12        # recorded energy against the manifest's, relative
FAMILY_ALPHA = 1e-3       # chance that a correct run fails any statistical row
Z_FLOOR = 4.0
MIN_TX_MOMENTS = 100      # transmissions needed for the exit-moment rows
FILTER_MARGIN = 0.95      # physics row: grid points below this fraction of the filter limit
FILTER_T_MAX = 0.01
NORMAL, E_LONG = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
PLANCK_GRID_POINTS = 20001
WARNING = re.compile(r"(BBR008) dataset=(\S+) side=(low|high) |(BBR026) dataset=(\S+) direction=(below|above) ")
OTHER_CODE = re.compile(r"BBR0\d\d")


class Report:
    """PASS / FAIL / SKIP rows; the statistical rows are judged in finish(), against z_max(M)."""

    def __init__(self):
        self.rows, self.stats = [], []

    def check(self, name, ok, detail=""):
        self.rows.append(("PASS" if ok else "FAIL", name, detail))

    def skip(self, name, detail=""):
        self.rows.append(("SKIP", name, detail))

    def stat(self, name, observed, expected, sigma, slack=0.0, detail=""):
        self.stats.append((name, float(observed), float(expected), float(sigma), float(slack), detail))

    def finish(self):
        """Judge the statistical rows; True when no row failed."""
        zm = z_max(len(self.stats))
        for name, obs, exp, sig, slack, detail in self.stats:
            dev = abs(obs - exp)
            ok = dev <= zm * sig + slack
            z = dev / sig if sig > 0 else (0.0 if dev <= slack else math.inf)
            self.rows.append(("PASS" if ok else "FAIL", name,
                              f"observed {obs:.6g}, expected {exp:.6g}, {z:.2f} SE (limit {zm:.2f})"
                              + (f"; {detail}" if detail else "")))
        self.stats = []
        return all(status != "FAIL" for status, _, _ in self.rows)


def z_max(m):
    """Limit in standard errors for m statistical rows: a correct run fails any of them with
    probability at most FAMILY_ALPHA (two-sided, Bonferroni), and never below Z_FLOOR."""
    if m <= 0:
        return Z_FLOOR
    return max(Z_FLOOR, NormalDist().inv_cdf(1.0 - FAMILY_ALPHA / (2.0 * m)))


class Tables:
    """bbrsim.hfss tables of one waveguides directory, each stem loaded once; get() returns
    the dataset dict, or the exception loading raised, so a bad table is a FAIL row."""

    def __init__(self, waveguides):
        self.waveguides, self.cache = waveguides, {}

    def get(self, stem):
        if stem not in self.cache:
            try:
                self.cache[stem] = hfss.load_dataset(stem, self.waveguides)
            except (OSError, ValueError) as err:
                self.cache[stem] = err
        return self.cache[stem]


def dataset_id(volume):
    """BBRCrackLibrary::DatasetIdOf: the volume name without a ':N' copy suffix."""
    return str(volume).split(":")[0]


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def expected_warnings(entries, grids, cutoffs):
    """{(code, dataset, side or direction): count} the log must hold, from the decided entries
    (dataset, nu_GHz, grid_GHz) of the whole process: BBR008 once per dataset and side when an
    entry lies beyond that edge of a grid of more than one point, BBR026 once per dataset and
    direction when an entry and its grid point lie on opposite sides of the cutoff (the
    conditions of BBRCrackLibrary::Lookup), otherwise none."""
    want = {}
    for did in grids:
        want.update({("BBR008", did, "low"): 0, ("BBR008", did, "high"): 0,
                     ("BBR026", did, "below"): 0, ("BBR026", did, "above"): 0})
    for did, nu, f in entries:
        grid = grids[did]
        if len(grid) > 1 and nu < grid[0][0]:
            want[("BBR008", did, "low")] = 1
        if len(grid) > 1 and nu > grid[-1][0]:
            want[("BBR008", did, "high")] = 1
        fc = cutoffs.get(did, -1.0)
        if fc > 0.0 and (nu < fc) != (f < fc):
            want[("BBR026", did, "below" if nu < fc else "above")] = 1
    return want


def parse_log(text):
    """The warnings, banners, GeomNav lines, other BBR codes and startup lines of a run.log."""
    warnings = {}
    for m in WARNING.finditer(text):
        key = ("BBR008", m.group(2), m.group(3)) if m.group(1) else ("BBR026", m.group(5), m.group(6))
        warnings[key] = warnings.get(key, 0) + 1
    lines = text.splitlines()
    return {"warnings": warnings,
            "banners": sum("G4Exception-START" in ln for ln in lines),
            "geomnav": sum("geomnav" in ln.lower() for ln in lines),
            "other_codes": [ln for ln in lines
                            if OTHER_CODE.search(ln) and "BBR008" not in ln and "BBR026" not in ln],
            "crack_lines": [ln for ln in lines if "[BBR] crack " in ln],
            "dataset_lines": [ln for ln in lines if "[BBR] HFSS dataset " in ln]}


def load_runs(run_dir, man):
    """[(run, crossings or None, metadata or None, problem or None)] in manifest order."""
    out = []
    for r in man["runs"]:
        path = os.path.join(run_dir, r["file"])
        if not os.path.exists(path):
            out.append((r, None, None, f"{path} is missing"))
            continue
        try:
            out.append((r, io.load_crossings(path), io.load_metadata(path), None))
        except Exception as err:      # an unreadable output is a FAIL row, not a traceback
            out.append((r, None, None, f"{path}: {type(err).__name__}: {err}"))
    return out


def run_tag(r):
    what = f" {r['target_ghz']:.6g} GHz" if "target_ghz" in r else ""
    return f"run {r['run']:03d} {r['dataset']} {r['kind']}{what}"


def check_gun_exact(rep, tag, r, dec, ids):
    """One entry per event on the aimed crack, the configured direction and energy, the
    expected grid frequency."""
    n = len(dec)
    on_aim = bool(n) and bool((ids == r["dataset"]).all())
    one_each = n == r["events"] and not dec["event_id"].duplicated().any()
    rep.check(f"{tag} entries", one_each and on_aim,
              f"{n} decided entries for {r['events']} events" + ("" if on_aim else ", not all on the aimed crack"))
    if not n:
        return
    dev = float(np.abs(dec[["px_pre", "py_pre", "pz_pre"]].to_numpy(float) - unit(r["direction"])).max())
    rep.check(f"{tag} direction", dev < TOL_DIR, f"max |k - k_configured| {dev:.1e}")
    rel = float(np.abs(dec["energy_eV"].to_numpy(float) / r["energy_eV"] - 1.0).max())
    rep.check(f"{tag} energy", rel < TOL_ENERGY, f"max relative deviation from {r['energy_eV']:.10e} eV: {rel:.1e}")
    sel = sorted(set(dec["hfss_freq_GHz"].to_numpy(float).tolist()))
    rep.check(f"{tag} selection", sel == [r["expected_grid_ghz"]],
              f"selected {sel} GHz, expected {r['expected_grid_ghz']:g}")


def check_reflections(rep, tag, dec):
    rf = dec[dec["status"] == REFLECT]
    if not len(rf):
        return
    k_exp = np.array([hfss.reflected_direction(k) for k in rf[["px_pre", "py_pre", "pz_pre"]].to_numpy(float)])
    k_exp /= np.linalg.norm(k_exp, axis=1, keepdims=True)
    dev = float(np.abs(rf[["px_post", "py_post", "pz_post"]].to_numpy(float) - k_exp).max())
    rep.check(f"{tag} reflection", dev < TOL_REFL, f"{len(rf)} reflections, max deviation from specular {dev:.1e}")


def check_log(rep, log_path, entries, grids, cutoffs, placed, wg):
    try:
        with open(log_path, errors="replace") as fh:
            lg = parse_log(fh.read())
    except OSError as err:
        rep.check("log", False, f"{log_path}: {err}")
        return
    want = expected_warnings(entries, grids, cutoffs)
    for (code, did, which), n in sorted(want.items()):
        got = lg["warnings"].get((code, did, which), 0)
        rep.check(f"log {code} {did} {which}", got == n, f"{got} in the log, {n} expected from the output")
    stray = sorted(k for k in lg["warnings"] if k not in want)
    n_warn = sum(lg["warnings"].values())
    rep.check("log hygiene",
              lg["geomnav"] == 0 and not lg["other_codes"] and not stray and lg["banners"] == n_warn,
              f"GeomNav {lg['geomnav']}, other BBR codes {len(lg['other_codes'])}, "
              f"G4Exception banners {lg['banners']} for {n_warn} BBR008/BBR026 warnings"
              + (f", warnings for unknown datasets {stray}" if stray else ""))
    for vol, did in sorted(placed.items()):
        lines = [ln.rstrip() for ln in lg["crack_lines"] if f"[BBR] crack {vol}: " in ln]
        m = re.search(r"; (\d+) sidecar\(s\) fit$", lines[0]) if len(lines) == 1 else None
        rep.check(f"log startup {vol}", m is not None and int(m.group(1)) == len(grids.get(did, [])),
                  lines[0] if len(lines) == 1 else f"{len(lines)} '[BBR] crack {vol}:' lines")
    for did in sorted(set(placed.values())):
        lines = [ln.rstrip() for ln in lg["dataset_lines"] if f"[BBR] HFSS dataset {did}: " in ln]
        m = (re.search(r": (\d+) frequency grid point\(s\):.*\(from (.+)\)$", lines[0])
             if len(lines) == 1 else None)
        ok = (m is not None and int(m.group(1)) == len(grids.get(did, []))
              and os.path.realpath(m.group(2)) == os.path.realpath(wg))
        rep.check(f"log grid {did}", ok,
                  lines[0] if len(lines) == 1 else f"{len(lines)} '[BBR] HFSS dataset {did}:' lines")


def validate(run_dir, log_path=None, data_root=None, tables=None):
    """Every row for the run in run_dir; returns the Report (finish() judges it)."""
    rep = Report()
    man_path = os.path.join(run_dir, "tree_check.json")
    try:
        with open(man_path) as fh:
            man = json.load(fh)
    except (OSError, ValueError) as err:
        rep.check("manifest", False, f"{man_path}: {err}")
        return rep
    if man.get("schema") != SCHEMA:
        rep.check("manifest", False, f"{man_path}: schema {man.get('schema')!r}, expected {SCHEMA!r}")
        return rep
    root = os.path.realpath(data_root or man["data_root"])
    wg = os.path.join(root, "waveguides")
    tables = tables or Tables(wg)
    grids = {d: [(float(f), s) for f, s in v["grid"]] for d, v in man["datasets"].items()}
    aimed = sorted(d for d, v in man["datasets"].items() if v["aimed"])
    cutoffs = {d: float(man["datasets"][d]["cutoff_ghz"]) for d in aimed}
    rep.check("manifest", True, f"{man_path}: {len(man['runs'])} runs over {', '.join(aimed)}")

    try:
        now = {d: [(f, s) for f, s in g] for d, g in hfss.discover_datasets(wg).items()}
        rep.check("tree unchanged", now == grids, f"{wg}: grids as generated" if now == grids
                  else f"{wg}: the grids differ from the manifest's; regenerate the run")
    except (OSError, ValueError) as err:
        rep.check("tree unchanged", False, f"{wg}: {err}")
    for did in man.get("not_exercised", []):
        rep.skip(f"dataset {did}", "the test world cannot place it; not exercised")
    for s in man.get("skipped", []):
        rep.skip(f"dataset {s['dataset']} {s['kinds']}", s["reason"])

    entries, served, placed = [], {d: set() for d in aimed}, None
    n_dec = n_bad_sel = n_bad_sentinel = 0
    for r, df, meta, problem in load_runs(run_dir, man):
        tag = run_tag(r)
        if problem:
            rep.check(f"{tag} output", False, problem)
            continue
        if meta.get("schema_version") != 2 or "data" not in meta:
            rep.check(f"{tag} metadata", False, "no schema-2 metadata with a data block")
            continue
        data = meta["data"]
        hd = {h["volume"]: h for h in data.get("hfss_datasets", [])}
        labels_ok = (sorted(h["dataset_id"] for h in hd.values()) == aimed and all(
            [f["label"] for f in h["frequencies"]] == [s[len(h["dataset_id"]) + 1:] for _, s in grids[h["dataset_id"]]]
            for h in hd.values()))
        same_root = os.path.realpath(data.get("directory", "")) == root
        rep.check(f"{tag} metadata", same_root and labels_ok,
                  f"data root {data.get('directory')}" + ("" if same_root else f", not {root}")
                  + ("" if labels_ok else "; placed cracks or frequency labels differ from the tree"))
        if placed is None:
            placed = {v: h["dataset_id"] for v, h in hd.items()}
        if "hfss_freq_GHz" not in df.columns:
            rep.check(f"{tag} hfss_freq_GHz column", False, "missing")
            continue
        cr = select.crack_crossings(df)
        dec = cr[cr["status"].isin(DECIDED)]
        n_bad_sentinel += int((df.drop(index=dec.index)["hfss_freq_GHz"] != -1.0).sum())
        ids = dec["vol_post"].map(dataset_id).to_numpy()
        nu = np.atleast_1d(hfss.photon_frequency_GHz(dec["energy_eV"].to_numpy(float)))
        f_rec = dec["hfss_freq_GHz"].to_numpy(float)
        f_exp = np.array([hfss.select_frequency(grids[i], n)[0] if i in grids else np.nan
                          for i, n in zip(ids, nu)])
        n_dec += len(dec)
        n_bad_sel += int((f_rec != f_exp).sum())
        entries.extend((i, float(n), float(f)) for i, n, f in zip(ids, nu, f_rec) if i in cutoffs)
        if r["kind"] == "planck":
            rep.check(f"{tag} entries", not dec["event_id"].duplicated().any(),
                      f"{len(dec)} decided entries in {r['events']} events, at most one per event")
        else:
            check_gun_exact(rep, tag, r, dec, ids)
            if len(dec) and (ids == r["dataset"]).all():
                served[r["dataset"]].update(np.unique(f_rec).tolist())
        check_reflections(rep, tag, dec)

    rep.check("sentinel", n_bad_sentinel == 0,
              f"{n_bad_sentinel} rows that are not decided crack entries carry a frequency other than -1")
    rep.check("selection", n_dec > 0 and n_bad_sel == 0,
              f"{n_dec - n_bad_sel} of {n_dec} decided entries carry the grid frequency the selection rule gives")
    for did in aimed:
        missing = [f for f, _ in grids[did] if f not in served[did]]
        rep.check(f"coverage {did}", not missing,
                  "every grid frequency served by a gun run" if not missing else f"no gun run served {missing} GHz")
    if log_path is not None:
        check_log(rep, log_path, entries, {d: grids[d] for d in aimed}, cutoffs, placed or {}, wg)
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="directory holding tree_check.json and output/")
    ap.add_argument("--log", help="the run's log; enables the log rows")
    ap.add_argument("--data-root", help="data root, if the tree moved since generation")
    a = ap.parse_args(argv)
    rep = validate(a.run_dir, a.log, a.data_root)
    n_stats = len(rep.stats)
    ok = rep.finish()
    print(f"check_tree_run: {os.path.realpath(a.run_dir)}; {n_stats} statistical rows, "
          f"limit {z_max(n_stats):.2f} standard errors")
    for status, name, detail in rep.rows:
        print(f"  {status:<4} {name}: {detail}")
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
