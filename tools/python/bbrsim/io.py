"""Load BBRsim ROOT output (output/bbr.root) into tidy pandas DataFrames.

The C++ writes categorical fields (status, event_type, volume, material) as
integer codes; each ROOT file has a matching `<stem>.metadata.json` containing
its legend and provenance. Legacy files use sibling `bbr_legend.json`. This module decodes them to readable strings so downstream analysis
matches the old CSV column semantics.
"""
import json

import uproot


def load_metadata(root_path):
    """Load this result's metadata; legacy files expose only their shared legend.

    A present but corrupt/newer metadata file is an error, never a reason to
    silently substitute a possibly unrelated legacy legend.
    """
    from pathlib import Path
    path = Path(root_path).with_suffix(".metadata.json")
    if path.exists():
        with path.open() as fh:
            metadata = json.load(fh)
        if metadata.get("schema_version") != 2:
            raise ValueError(f"Unsupported BBRsim schema version in {path}")
        if not isinstance(metadata.get("legend"), dict):
            raise ValueError(f"Missing legend in {path}")
        return metadata
    legend_path = Path(root_path).parent / "bbr_legend.json"
    with legend_path.open() as fh:
        return {"schema_version": 1, "legend": json.load(fh)}


def _legend_for(root_path):
    raw = load_metadata(root_path)["legend"]
    return {cat: {int(k): v for k, v in m.items()} for cat, m in raw.items()}


def load(path):
    """Return (crossings_df, abspoints_df) with decoded string columns."""
    metadata = load_metadata(path)
    maps = {cat: {int(k): v for k, v in m.items()} for cat, m in metadata["legend"].items()}
    with uproot.open(path) as f:
        cr = f["crossings"].arrays(library="pd")
        ab = f["abspoints"].arrays(library="pd")

    if metadata["schema_version"] == 1 and "n_boundary" in cr.columns:
        raise ValueError("Schema 2 output is missing its per-result .metadata.json file")

    cr["status"] = cr["status_code"].map(maps["status"])
    cr["event_type"] = cr["event_type_code"].map(maps["event_type"])
    cr["vol_pre"] = cr["vol_pre_code"].map(maps["volume"])
    cr["mat_pre"] = cr["mat_pre_code"].map(maps["material"])
    cr["vol_post"] = cr["vol_post_code"].map(maps["volume"])
    cr["mat_post"] = cr["mat_post_code"].map(maps["material"])

    # Always present, empty when there are no rows (an empty tree keeps its columns).
    ab["term_vol"] = ab["term_vol_code"].map(maps["volume"])
    ab["term_status"] = ab["term_status_code"].map(maps["status"])
    return cr, ab


def load_crossings(path):
    """Convenience: crossings DataFrame with the legacy CSV column names."""
    cr, _ = load(path)
    return cr


def load_many(paths):
    """Concatenate the crossings of several ROOT files (one per run).

    Multi-run macros name each run's file with /analysis/setFileName because
    the run action reopens the output file at every /run/beamOn; the run_id
    column tells the runs apart after concatenation.
    """
    import pandas as pd
    frames = [load_crossings(p) for p in paths]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
