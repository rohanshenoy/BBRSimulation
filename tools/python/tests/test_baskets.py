"""bbrsim.io.basket_problems: the basket-integrity check of a ROOT file (G4-026, G4-035)."""
import uproot
from bbrsim import io

import basket_fixture as bf


def test_clean_files_pass(tmp_path):
    for i, compression in enumerate((None, uproot.ZLIB(1))):
        p = tmp_path / f"clean{i}.root"
        bf.write_tree(p, compression=compression)
        r = io.basket_problems(p)
        assert (r.problems, r.on_disk, r.embedded) == ([], 4, 0)    # two branches of two baskets each


def test_g4026_signature_is_reported(tmp_path):
    p = tmp_path / "g4026.root"
    bf.write_tree(p)
    bf.corrupt_g4026(p, "abspoints", "event_id", 0)
    with uproot.open(p) as f:
        ids = f["abspoints"]["event_id"].array(library="np")
    # The reader takes the block as data: "ZL\x08" and the low byte of c = 391, G4-026's 1514932415-style garbage.
    assert ids[0] == 0x5A4C0887 and list(ids[100:103]) == [100, 101, 102]
    r = io.basket_problems(p)
    assert r.on_disk == 4 and len(r.problems) == 1, r.problems
    assert r.problems[0].startswith("abspoints/event_id basket 0 (entries [0, 100), seek ")
    assert "keyed uncompressed but holds one whole ZL compression block of exactly its length" in r.problems[0]


def test_raw_data_that_starts_like_a_block_passes(tmp_path):
    """Raw data whose first bytes read "ZL\\x08" is not G4-026's basket unless the header's sizes fit the basket and,
    for ZL, the stream decompresses to them."""
    p = tmp_path / "zl_raw.root"
    bf.write_tree(p)
    bf.rewrite_payload(p, "abspoints", "event_id", 0, bf.zl_header(100, 400) + bytes(391))           # sizes disagree
    bf.rewrite_payload(p, "abspoints", "event_id", 1, bf.zl_header(391, 400) + bytes(range(256)) + bytes(135))
    assert io.basket_problems(p).problems == []                                    # sizes fit, not a zlib stream


def test_unreadable_basket_is_reported(tmp_path):
    p = tmp_path / "unreadable.root"
    bf.write_tree(p, compression=uproot.ZLIB(1))
    _seek, nbytes, _objlen, keylen = bf.basket_key(p, "abspoints", "x_mm", 0)
    bf.rewrite_payload(p, "abspoints", "x_mm", 0, b"ZL\x08" + bytes(nbytes - keylen - 3))   # a block header over zeros
    problems = io.basket_problems(p).problems
    assert len(problems) == 1 and problems[0].startswith("abspoints/x_mm cannot be read: "), problems


class FakeBranch:
    """The members basket_problems reads from a TBranch."""
    def __init__(self, members, num_baskets):
        self.members, self.num_baskets, self.name = members, num_baskets, "event_id"

    def member(self, key):
        return self.members[key]


def test_embedded_baskets_are_counted_not_failed(tmp_path):
    """A basket kept in the TBranch object has no on-disk key: it is skipped for the signature (the decode check
    covers it) and counted; it never fails for want of a seek."""
    p = tmp_path / "embedded.root"
    bf.write_tree(p, baskets=1)
    seek = bf.basket_key(p, "abspoints", "event_id", 0)[0]
    with open(p, "rb") as fh:
        size = fh.seek(0, 2)
        # one written basket at seek 0 (no key), one at its real seek, one more embedded past fWriteBasket
        fake = FakeBranch({"fWriteBasket": 2, "fBasketSeek": [0, seek, 0], "fBasketEntry": [0, 100, 200, 300]}, 3)
        problems, on_disk, embedded = io._branch_basket_problems(fh, size, "abspoints", fake)
    assert (problems, on_disk, embedded) == ([], 1, 2)


def test_unreadable_basket_list_fails(tmp_path):
    p = tmp_path / "nolist.root"
    bf.write_tree(p, baskets=1)
    with open(p, "rb") as fh:
        size = fh.seek(0, 2)
        problems, on_disk, embedded = io._branch_basket_problems(fh, size, "abspoints", FakeBranch({}, 1))
    assert on_disk == 0 and len(problems) == 1
    assert problems[0].startswith("abspoints/event_id: fWriteBasket or fBasketSeek cannot be read")


def test_key_past_the_end_of_the_file_fails(tmp_path):
    p = tmp_path / "short.root"
    bf.write_tree(p, baskets=1)
    with open(p, "rb") as fh:
        size = fh.seek(0, 2)
        fake = FakeBranch({"fWriteBasket": 1, "fBasketSeek": [size + 10], "fBasketEntry": [0, 100]}, 1)
        problems, _on_disk, _embedded = io._branch_basket_problems(fh, size, "abspoints", fake)
    assert len(problems) == 1 and "past the end of the file" in problems[0], problems
