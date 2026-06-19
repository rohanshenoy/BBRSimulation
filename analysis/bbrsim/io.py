"""Load BBRsim ROOT output (output/bbr.root) into tidy pandas DataFrames.

The C++ writes categorical fields (status, event_type, volume, material) as
integer codes; a sibling `bbr_legend.json` ({category: {code: name}}) maps them
back. This module decodes them to readable strings so downstream analysis
matches the old CSV column semantics.
"""
import json
import os

import uproot


def _legend_for(root_path):
    """Load the bbr_legend.json sidecar that sits next to the ROOT file."""
    legend_path = os.path.join(os.path.dirname(root_path), "bbr_legend.json")
    with open(legend_path) as fh:
        raw = json.load(fh)
    # JSON keys are strings; convert code keys to int.
    return {cat: {int(k): v for k, v in m.items()} for cat, m in raw.items()}


def load(path):
    """Return (crossings_df, abspoints_df) with decoded string columns."""
    f = uproot.open(path)
    maps = _legend_for(path)

    cr = f["crossings"].arrays(library="pd")
    cr["status"] = cr["status_code"].map(maps["status"])
    cr["event_type"] = cr["event_type_code"].map(maps["event_type"])
    cr["vol_pre"] = cr["vol_pre_code"].map(maps["volume"])
    cr["mat_pre"] = cr["mat_pre_code"].map(maps["material"])
    cr["vol_post"] = cr["vol_post_code"].map(maps["volume"])
    cr["mat_post"] = cr["mat_post_code"].map(maps["material"])

    ab = f["abspoints"].arrays(library="pd")
    if len(ab):
        ab["term_vol"] = ab["term_vol_code"].map(maps["volume"])
        ab["term_status"] = ab["term_status_code"].map(maps["status"])
    return cr, ab


def load_crossings(path):
    """Convenience: crossings DataFrame with the legacy CSV column names."""
    cr, _ = load(path)
    return cr
