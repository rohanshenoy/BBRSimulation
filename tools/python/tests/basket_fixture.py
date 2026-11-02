"""Synthetic ROOT files for the basket-integrity tests (G4-026, G4-035).

G4-026: g4tools before Geant4 11.3.0 could compress a basket to exactly its raw length and then key it as
uncompressed (ObjLen == Nbytes - KeyLen), so every reader takes the compression block for the data. These helpers
write small TTrees with uproot and rewrite one basket's payload in place to reproduce that, or other faults.
"""
import struct
import zlib

import numpy as np
import uproot


def write_tree(path, rows=100, baskets=2, compression=None, name="abspoints"):
    """A TTree of event_id (int32) and x_mm (float64), one basket per extend of `rows` rows."""
    with uproot.recreate(path, compression=compression) as f:
        tree = f.mktree(name, {"event_id": np.int32, "x_mm": np.float64})
        for b in range(baskets):
            tree.extend({"event_id": np.arange(b * rows, (b + 1) * rows, dtype=np.int32), "x_mm": np.ones(rows)})


def basket_key(path, tree, branch, i):
    """(seek, Nbytes, ObjLen, KeyLen) of on-disk basket i of tree/branch."""
    with uproot.open(path) as f:
        seek = int(f[tree][branch].member("fBasketSeek")[i])
    with open(path, "rb") as fh:
        fh.seek(seek)
        nbytes, _version, objlen, _datime, keylen, _cycle = struct.unpack(">ihiIhh", fh.read(18))
    return seek, nbytes, objlen, keylen


def rewrite_payload(path, tree, branch, i, payload):
    """Overwrite basket i's payload in place; the payload must keep the basket's on-disk length."""
    seek, nbytes, _objlen, keylen = basket_key(path, tree, branch, i)
    assert len(payload) == nbytes - keylen, (len(payload), nbytes - keylen)
    with open(path, "r+b") as fh:
        fh.seek(seek + keylen)
        fh.write(payload)


def zl_header(c, u):
    """ROOT's 9-byte zlib block header: "ZL", method 8, compressed and uncompressed sizes (3 bytes each, little-endian)."""
    return b"ZL\x08" + c.to_bytes(3, "little") + u.to_bytes(3, "little")


def zl_block(length):
    """(raw, block): a ROOT zlib block of exactly `length` bytes that decompresses to `raw`, also `length` bytes long:
    the compressed-to-its-own-length case of G4-026. raw is int32 zeros after a run of k random values; the
    compressed length grows with k almost monotonically, so k is found by bisection and the neighbours scanned."""
    assert length % 4 == 0 and length >= 64, length
    n, target = length // 4, length - 9
    for seed in range(50):
        values = np.random.default_rng(seed).integers(0, 2**31 - 1, n)

        def raw(k):
            v = np.zeros(n, dtype=">i4")
            v[:k] = values[:k]
            return v.tobytes()

        for level in (1, 6, 9):
            lo, hi = 0, n
            while lo < hi:
                mid = (lo + hi) // 2
                if len(zlib.compress(raw(mid), level)) < target:
                    lo = mid + 1
                else:
                    hi = mid
            for k in range(max(0, lo - 40), min(n, lo + 40) + 1):
                r = raw(k)
                z = zlib.compress(r, level)
                if len(z) == target:
                    return r, zl_header(len(z), length) + z
    raise RuntimeError(f"no zlib stream of {target} bytes found")


def corrupt_g4026(path, tree, branch, i=0):
    """Rewrite uncompressed basket i of tree/branch as G4-026 left it: a whole zlib block of exactly the raw length
    under a key that still says uncompressed. Returns the raw bytes the block decompresses to."""
    _seek, nbytes, objlen, keylen = basket_key(path, tree, branch, i)
    assert objlen == nbytes - keylen, "the basket must be keyed uncompressed"
    raw, block = zl_block(objlen)
    rewrite_payload(path, tree, branch, i, block)
    return raw
