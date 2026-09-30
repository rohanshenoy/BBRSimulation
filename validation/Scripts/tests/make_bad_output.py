"""
make_bad_output.py
Write a small synthetic BBRsim output (bbr.root + bbr_legend.json) for the
validator negative tests: a clean file that check_invariants.py passes, and
switches that each break one invariant or the legend.

    conda run -n bbrsim python validation/Scripts/tests/make_bad_output.py OUT_DIR [switches]

    (no switch)                     clean: check_invariants.py PASS
    --metal-start                   one crossings row starts in Cu (mat_pre = Cu_RRR100_T4K)
    --drop-legend-code CAT:CODE     remove CODE from legend category CAT (e.g. material:0,
                                    the Cu code; a code the data use then decodes to NaN)
    --unknown-term                  one more abspoints row, labelled "unknown"
    --empty                         no crossings rows (only world exits in abspoints, so
                                    --allow-no-crossings makes it pass)

The trees are written with uproot mktree + extend (a dict assignment would write
an RNTuple, which bbrsim.io cannot read). Nothing is written outside OUT_DIR.
"""
import argparse
import json
import os
import sys

import numpy as np
import uproot

# The code maps of a real test-world run (output/bbr_legend.json).
LEGEND = {
    "status": {"0": "FresnelRefraction", "1": "FresnelReflection", "2": "TIR",
               "3": "LambertianReflection", "4": "LobeReflection", "5": "SpikeReflection",
               "6": "BackScattering", "7": "Absorption", "8": "Detection", "9": "NotAtBoundary",
               "10": "SameMaterial", "11": "StepTooSmall", "12": "NoRINDEX", "13": "Other",
               "14": "BBRDiffractionTransmit", "15": "BBRDiffractionReflect", "16": "BBRReflect",
               "17": "BBRAbsorb", "18": "unknown", "19": "WorldExit", "20": "BulkAbsorption"},
    "event_type": {"0": "transmission", "1": "reflection", "2": "absorption", "3": "other"},
    "volume": {"0": "CuSlab", "1": "InfParallelPlate_crack1Rohan", "2": "InfParallelPlate_crack2",
               "3": "World", "-1": "none"},
    "material": {"0": "Cu_RRR100_T4K", "1": "G4_Galactic", "2": "vacuum_wg", "-1": "none"},
}
CU, GALACTIC, VACUUM_WG = 0, 1, 2          # material codes
CUSLAB, CRACK1, WORLD, NONE = 0, 1, 3, -1  # volume codes
REFLECT, ABSORB, TRANSMIT, UNKNOWN, WORLDEXIT = 16, 17, 14, 18, 19   # status codes

I4, F8 = np.int32, np.float64
CR_INT = ["run_id", "event_id", "vol_pre_code", "mat_pre_code", "vol_post_code", "mat_post_code",
          "status_code", "event_type_code", "n_reflect"]
CR_DBL = ["x_mm", "y_mm", "z_mm", "energy_eV", "px_pre", "py_pre", "pz_pre", "px_post", "py_post",
          "pz_post", "theta_in_deg", "phi_in_deg", "hfss_freq_GHz"]
AB_INT = ["run_id", "event_id", "n_reflect", "term_vol_code", "term_status_code"]
AB_DBL = ["x_mm", "y_mm", "z_mm", "energy_eV", "px", "py", "pz"]

# Clean content, World-side as the stepping action logs it. Event 0 reflects off
# the Cu and is then absorbed by it; event 1 enters crack1 and leaves the world;
# event 2 leaves the world without crossing anything.
CROSSINGS = {"event_id":        [0, 0, 1],
             "n_reflect":       [1, 2, 1],
             "vol_pre_code":    [WORLD, WORLD, WORLD],
             "mat_pre_code":    [GALACTIC, GALACTIC, GALACTIC],
             "vol_post_code":   [CUSLAB, CUSLAB, CRACK1],
             "mat_post_code":   [CU, CU, VACUUM_WG],
             "status_code":     [REFLECT, ABSORB, TRANSMIT],
             "event_type_code": [1, 2, 0],
             "hfss_freq_GHz":   [-1., -1., 500.]}
ABSPOINTS = {"event_id":         [0, 1, 2],
             "n_reflect":        [2, 1, 0],
             "term_vol_code":    [CUSLAB, NONE, NONE],
             "term_status_code": [ABSORB, WORLDEXIT, WORLDEXIT]}


def columns(content, ints, dbls):
    n = len(next(iter(content.values())))
    out = {c: np.asarray(content.get(c, [0] * n), I4) for c in ints}
    out.update({c: np.asarray(content.get(c, [0.] * n), F8) for c in dbls})
    out["energy_eV"] = np.full(n, 2.07e-3, F8)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out_dir")
    ap.add_argument("--metal-start", action="store_true")
    ap.add_argument("--drop-legend-code", action="append", default=[], metavar="CAT:CODE")
    ap.add_argument("--unknown-term", action="store_true")
    ap.add_argument("--empty", action="store_true")
    a = ap.parse_args()

    cr = {k: list(v) for k, v in CROSSINGS.items()}
    ab = {k: list(v) for k, v in ABSPOINTS.items()}
    if a.metal_start:
        cr["vol_pre_code"][0], cr["mat_pre_code"][0] = CUSLAB, CU
    if a.unknown_term:
        for k, v in (("event_id", 3), ("n_reflect", 1), ("term_vol_code", CUSLAB),
                     ("term_status_code", UNKNOWN)):
            ab[k].append(v)
    if a.empty:
        cr = {k: [] for k in cr}
        keep = [i for i, s in enumerate(ab["term_status_code"]) if s == WORLDEXIT]
        ab = {k: [v[i] for i in keep] for k, v in ab.items()}
    legend = json.loads(json.dumps(LEGEND))
    for spec in a.drop_legend_code:
        cat, _, code = spec.partition(":")
        if cat not in legend or code not in legend[cat]:
            ap.error(f"--drop-legend-code {spec}: no code {code!r} in category {cat!r}")
        del legend[cat][code]

    crossings = columns(cr, CR_INT, CR_DBL)
    abspoints = columns(ab, AB_INT, AB_DBL)
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, "bbr.root")
    with uproot.recreate(path) as f:
        for name, arr in (("crossings", crossings), ("abspoints", abspoints)):
            tree = f.mktree(name, {k: v.dtype for k, v in arr.items()})
            if len(arr["event_id"]):
                tree.extend(arr)
    with open(os.path.join(a.out_dir, "bbr_legend.json"), "w") as fh:
        json.dump(legend, fh)
    print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
