"""Validate per-track boundary and reflection counters in BBRsim ROOT output.

Both counters include the current crossing. A world exit without any logged
boundary has zero for each. The legacy n_reflect column counts crossings, not
reflections, and is deliberately ignored here.

Usage: conda run -n bbrsim python validation/check_nreflect.py [path/to/bbr.root]
"""

import argparse
import sys

import numpy as np

from bbrsim.io import load


REFLECTIONS = {
    "FresnelReflection", "TIR", "LambertianReflection", "LobeReflection",
    "SpikeReflection", "BackScattering", "BBRDiffractionReflect", "BBRReflect",
    "CoatedDielectricReflection",
    "PolishedLumirrorAirReflection", "PolishedLumirrorGlueReflection",
    "PolishedAirReflection", "PolishedTeflonAirReflection",
    "PolishedTiOAirReflection", "PolishedTyvekAirReflection",
    "PolishedVM2000AirReflection", "PolishedVM2000GlueReflection",
    "EtchedLumirrorAirReflection", "EtchedLumirrorGlueReflection",
    "EtchedAirReflection", "EtchedTeflonAirReflection",
    "EtchedTiOAirReflection", "EtchedTyvekAirReflection",
    "EtchedVM2000AirReflection", "EtchedVM2000GlueReflection",
    "GroundLumirrorAirReflection", "GroundLumirrorGlueReflection",
    "GroundAirReflection", "GroundTeflonAirReflection",
    "GroundTiOAirReflection", "GroundTyvekAirReflection",
    "GroundVM2000AirReflection", "GroundVM2000GlueReflection",
}
OTHER_BOUNDARY_STATUSES = {
    "FresnelRefraction", "Absorption", "Detection", "NotAtBoundary",
    "SameMaterial", "StepTooSmall", "NoRINDEX",
    "BBRDiffractionTransmit", "BBRAbsorb",
    "Transmission", "Dichroic", "CoatedDielectricRefraction",
    "CoatedDielectricFrustratedTransmission",
}
VALID_BOUNDARY_STATUSES = REFLECTIONS | OTHER_BOUNDARY_STATUSES
KILLING_BOUNDARY_STATUSES = {"BBRAbsorb", "Absorption", "Detection", "NoRINDEX"}
KEY = ("run_id", "event_id", "track_id")
REQUIRED = (*KEY, "n_boundary", "n_reflections")


def validate(crossings, abspoints):
    """Return (final reflection counts, errors) for all terminated tracks."""
    missing = {
        name: sorted(set(REQUIRED) - set(frame.columns))
        for name, frame in (("crossings", crossings), ("abspoints", abspoints))
    }
    if any(missing.values()):
        return [], ["unsupported schema: missing " + "; ".join(
            f"{name} {columns}" for name, columns in missing.items() if columns)]
    if abspoints.empty:
        return [], ["empty output: no abspoints termination rows"]

    errors = []
    counts = []
    term_groups = {}
    for key, group in abspoints.groupby(list(KEY), sort=False):
        term_groups[key] = group
        if len(group) != 1:
            errors.append(f"{key}: duplicate termination rows ({len(group)})")

    cross_groups = {key: group for key, group in crossings.groupby(list(KEY), sort=False)}
    for key in cross_groups.keys() - term_groups.keys():
        errors.append(f"{key}: missing termination row")

    for key, term_group in term_groups.items():
        if len(term_group) != 1:
            continue
        term = term_group.iloc[0]
        if not isinstance(term["term_status"], str) or term["term_status"] == "unknown":
            errors.append(f"{key}: unknown termination status {term['term_status']!r}")

        group = cross_groups.get(key)
        if group is None:
            if term["term_status"] in KILLING_BOUNDARY_STATUSES:
                errors.append(f"{key}: missing killing boundary crossing")
            if int(term["n_boundary"]) != 0 or int(term["n_reflections"]) != 0:
                errors.append(f"{key}: termination totals must be zero without crossings")
            counts.append(int(term["n_reflections"]))
            continue

        group = group.sort_values("n_boundary")
        if term["term_status"] in KILLING_BOUNDARY_STATUSES and group.iloc[-1]["status"] != term["term_status"]:
            errors.append(f"{key}: killing boundary status {term['term_status']!r} "
                          f"does not match last crossing {group.iloc[-1]['status']!r}")
        ordinals = group["n_boundary"].tolist()
        if len(ordinals) != len(set(ordinals)):
            errors.append(f"{key}: duplicate boundary ordinal")
        if ordinals != list(range(1, len(group) + 1)):
            errors.append(f"{key}: gap or invalid boundary ordinal sequence {ordinals}")

        reflected = 0
        for row in group.itertuples(index=False):
            if row.status not in VALID_BOUNDARY_STATUSES:
                errors.append(f"{key}: unknown status {row.status!r} at boundary {row.n_boundary}")
                continue
            reflected += row.status in REFLECTIONS
            if int(row.n_reflections) != reflected:
                errors.append(f"{key}: reflection count at boundary {row.n_boundary}: "
                              f"recorded {row.n_reflections}, expected {reflected}")

        if int(term["n_boundary"]) != len(group) or int(term["n_reflections"]) != reflected:
            errors.append(f"{key}: termination totals ({term['n_boundary']}, "
                          f"{term['n_reflections']}) disagree with crossings "
                          f"({len(group)}, {reflected})")
        counts.append(int(term["n_reflections"]))
    return counts, errors


def plot_counts(counts, output):
    """Plot the number of terminated tracks at each actual reflection count."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values, frequencies = np.unique(counts, return_counts=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(values, frequencies, width=0.85, color="steelblue")
    ax.set_xlabel("Actual reflections per track")
    ax.set_ylabel("Terminated tracks")
    ax.set_title(f"Per-track reflection count (N_tracks={len(counts)})")
    ax.set_xticks(values)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", default="output/bbr.root")
    args = parser.parse_args()

    crossings, abspoints = load(args.path)
    counts, errors = validate(crossings, abspoints)
    print(f"crossing rows : {len(crossings)}")
    print(f"terminated tracks : {len(abspoints)}")
    for error in errors[:20]:
        print(f"  {error}")
    if len(errors) > 20:
        print(f"  ... and {len(errors) - 20} more errors")
    if errors:
        print("RESULT: FAIL")
        return 1

    print(f"reflection counts : {dict(zip(*np.unique(counts, return_counts=True)))}")
    output = "nreflect_distribution.png"
    plot_counts(counts, output)
    print(f"Plot saved: {output}")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
