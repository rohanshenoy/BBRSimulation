"""Load BBRsim ROOT output (output/bbr.root) into tidy pandas DataFrames.

The C++ writes categorical fields (status, event_type, volume, material) as
integer codes; each ROOT file has a matching `<stem>.metadata.json` containing
its legend and provenance. Legacy files use sibling `bbr_legend.json`. This module decodes them to readable strings so downstream analysis
matches the old CSV column semantics.
"""
import json
import struct
import zlib
from dataclasses import dataclass

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


# --- Basket integrity (G4-026, G4-035) -------------------------------------------------------------------------

# ROOT compression block magics: zlib, LZMA, LZ4, ZSTD, and the old ROOT algorithm.
_BLOCK_MAGICS = (b"ZL", b"XZ", b"L4", b"ZS", b"CS")


@dataclass
class BasketCheck:
    """What basket_problems found: one line per problem; the baskets read through their on-disk keys; and the
    baskets kept in a TBranch object (no key on disk), which only the decode check covers."""
    problems: list
    on_disk: int
    embedded: int


def _whole_block(payload, objlen):
    """The magic of a payload that is one whole ROOT compression block decompressing to objlen bytes, or None.
    The 9-byte header gives the compressed size c and the uncompressed size u (3 bytes each, little-endian);
    c + 9 must be the payload's length and u the basket's ObjLen, and a ZL block must zlib-decompress to u
    bytes, so raw data that merely starts with a magic does not count."""
    magic = bytes(payload[:2])
    if len(payload) < 9 or magic not in _BLOCK_MAGICS:
        return None
    c, u = int.from_bytes(payload[3:6], "little"), int.from_bytes(payload[6:9], "little")
    if c + 9 != len(payload) or u != objlen:
        return None
    if magic == b"ZL":
        try:
            if len(zlib.decompress(payload[9:])) != u:
                return None
        except zlib.error:
            return None
    return magic.decode()


def _branch_basket_problems(fh, size, tree_name, branch):
    """(problems, on_disk, embedded) for one branch: the signature check of basket_problems on every basket
    with an on-disk key (fh, of `size` bytes, is the open file)."""
    where = f"{tree_name}/{branch.name}"
    try:
        n = int(branch.member("fWriteBasket"))
        seeks = [int(s) for s in branch.member("fBasketSeek")[:n]]
    except Exception as err:
        return [f"{where}: fWriteBasket or fBasketSeek cannot be read ({type(err).__name__}: {err})"], 0, 0
    try:
        entries = [int(e) for e in branch.member("fBasketEntry")]
    except Exception:
        entries = []
    problems, on_disk = [], 0
    embedded = max(int(branch.num_baskets) - n, 0) + seeks.count(0)
    for i, seek in enumerate(seeks):
        if seek == 0:                       # kept in the TBranch object: no key to read
            continue
        on_disk += 1
        span = f"entries [{entries[i]}, {entries[i + 1]})" if i + 1 < len(entries) else "entries ?"
        at = f"{where} basket {i} ({span}, seek {seek})"
        if seek + 18 > size:
            problems.append(f"{at}: its key reaches past the end of the file ({size} bytes)")
            continue
        fh.seek(seek)
        nbytes, _version, objlen, _datime, keylen, _cycle = struct.unpack(">ihiIhh", fh.read(18))
        if seek + nbytes > size:
            problems.append(f"{at}: its key reaches past the end of the file ({size} bytes)")
            continue
        if objlen != nbytes - keylen:       # keyed compressed: the decode check covers it
            continue
        fh.seek(seek + keylen)
        magic = _whole_block(fh.read(objlen), objlen)
        if magic:
            problems.append(f"{at}: keyed uncompressed but holds one whole {magic} compression block of exactly its "
                            f"length (G4-026's signature), so readers take the compressed bytes for the data")
    return problems, on_disk, embedded


def basket_problems(path):
    """Check every basket of every TTree in a ROOT file; returns a BasketCheck.

    Two kinds of problem, one line each, naming the tree, the branch and the basket:
    - G4-026's signature (_whole_block): a basket keyed uncompressed (ObjLen == Nbytes - KeyLen) whose payload is
      one whole ROOT compression block of exactly its length. g4tools before Geant4 11.3.0 wrote such baskets
      when compression saved nothing, and every reader then takes the block for the data (an event_id of
      1514932415 is "ZL\\x08" and one more byte).
    - a branch uproot cannot read, or whose fWriteBasket or fBasketSeek cannot be read.
    Baskets kept in a TBranch object have no on-disk key: they are skipped by the signature check, counted, and
    covered by the read. The check is exact (no tolerance) and independent of the Geant4 version that wrote the
    file; its cost is one full read of every branch plus one key-header read per basket.
    """
    problems, on_disk, embedded = [], 0, 0
    with uproot.open(path) as f, open(path, "rb") as fh:
        size = fh.seek(0, 2)
        for tree_name, classname in f.classnames(cycle=False).items():
            if classname != "TTree":
                continue
            for branch in f[tree_name].itervalues(recursive=True):
                p, d, e = _branch_basket_problems(fh, size, tree_name, branch)
                problems += p
                on_disk += d
                embedded += e
                try:
                    branch.array(library="np")
                except Exception as err:
                    problems.append(f"{tree_name}/{branch.name} cannot be read: {type(err).__name__}: {err}")
    return BasketCheck(problems, on_disk, embedded)
