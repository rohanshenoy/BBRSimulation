"""Shared selections / reductions on decoded BBRsim crossings DataFrames.

Small helpers reused across the scripts/ validators and plots, built on the
string-decoded columns from bbrsim.io.load / load_crossings.
"""
import re

_CU_RE = re.compile(r"Cu_RRR(\d+)_T(\d+(?:\.\d+)?)K")


def add_evt_key(df):
    """Return df with an 'evt_key' column uniquely identifying a photon.

    Keys on (run_id, event_id) when run_id is present (multi-run sessions),
    else event_id alone. Mutates and returns df.
    """
    if "run_id" in df.columns:
        df["evt_key"] = list(zip(df["run_id"], df["event_id"]))
    else:
        df["evt_key"] = df["event_id"]
    return df


def cu_boundary(df):
    """Rows where Cu is the post-step material (Drude name 'Cu_RRR{N}_T{T}K')."""
    return df[df["mat_post"].str.startswith("Cu_RRR", na=False)]


def first_hit_cu(df):
    """First-contact Cu rows: n_reflect == 1 and Cu is the post-step material."""
    return df[(df["n_reflect"] == 1)
              & df["mat_post"].str.startswith("Cu_RRR", na=False)]


def crack_crossings(df):
    """Rows where the pre-step material is the crack flag material 'vacuum_wg'."""
    return df[df["mat_pre"] == "vacuum_wg"]


def cu_absorption_stats(df):
    """(n_hit, n_absorbed) per-photon at the Cu wall.

    n_hit      = unique photons with a Cu boundary crossing.
    n_absorbed = unique photons whose Cu crossing has status 'BBRAbsorb'.
    Requires the BBR wrapper statuses (BBRAbsorb / BBRReflect) in the data.
    """
    d = add_evt_key(df.copy())
    cu = cu_boundary(d)
    n_hit = cu["evt_key"].nunique()
    n_absorbed = cu[cu["status"] == "BBRAbsorb"]["evt_key"].nunique()
    return n_hit, n_absorbed


def parse_cu_rrr_t(material_name):
    """Parse 'Cu_RRR100_T4K' -> (100, 4.0). Returns None if it doesn't match."""
    m = _CU_RE.search(str(material_name))
    if not m:
        return None
    return int(m.group(1)), float(m.group(2))
