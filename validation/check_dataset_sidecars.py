"""
check_dataset_sidecars.py
Every HFSS dataset under the given waveguides directories is complete and has a
schema-1 sidecar that passes bbrsim.sidecar.check_full. A dataset <stem> =
<id>_<freq>GHz is the trio <stem>_Ephi=0/, <stem>_Ephi=1/ and
<stem>.dataset.json; any other file in the directory (a tree manifest,
SHA256SUMS) is ignored. Per dataset:
  * all three parts are present (a lone part is a FAIL, so a lost Ephi=0
    directory cannot silently drop a grid frequency);
  * the checks BBRsim's C++ runs (F1-F11, F13; codes BBR024 / BBR025);
  * the full modes block re-derived from the cross-section;
  * the sha256, size and row count of the four CSVs;
  * the CSV checks C1-C5.
The frequencies of one ID in one directory must also agree on the
frequency-independent blocks; a FAIL names the blocks that differ.

Output: one PASS or FAIL line per dataset, a FAIL line for a directory argument
that cannot be read or holds no dataset, a count, then RESULT: PASS or
RESULT: FAIL. Exit code 0 on PASS, 1 on FAIL.

Usage:
    conda run -n bbrsim python validation/check_dataset_sidecars.py <waveguides-dir> [more ...]
"""
import os
import sys

from bbrsim import hfss, sidecar

PARTS = (("_Ephi=0", "dir"), ("_Ephi=1", "dir"), (".dataset.json", "file"))


def datasets(base):
    """{stem: set of the parts present}; only the three dataset parts count."""
    found = {}
    for name in os.listdir(base):
        full = os.path.join(base, name)
        for suffix, kind in PARTS:
            if name.endswith(suffix) and len(name) > len(suffix) and \
                    (os.path.isdir(full) if kind == "dir" else os.path.isfile(full)):
                found.setdefault(name[:-len(suffix)], set()).add(suffix)
    return found


def check_one(base, stem, present):
    """The checked sidecar, or raise with the reason."""
    missing = [stem + s + ("/" if kind == "dir" else "") for s, kind in PARTS if s not in present]
    if missing:
        raise ValueError(f"incomplete dataset: missing {', '.join(missing)}")
    if "_" not in stem:
        raise ValueError(f"stem {stem!r} is not <id>_<freq>GHz")
    dataset_id, label = stem.rsplit("_", 1)
    return sidecar.check_full(base, dataset_id, stem, hfss.parse_frequency_GHz(label) or -1.0)


dirs = sys.argv[1:] or [hfss.default_base_dir()]
n_sets, n_failing, n_bad_dirs = 0, 0, 0
for base in dirs:
    try:
        found = datasets(base)
    except OSError as e:
        n_bad_dirs += 1
        print(f"  FAIL {base}: cannot list the directory ({type(e).__name__}: {e})")
        continue
    if not found:
        n_bad_dirs += 1
        print(f"  FAIL {base}: no dataset found")
        continue
    status, by_id = {}, {}          # stem -> None (pass) or reason; dataset id -> [(stem, sidecar)]
    for stem in sorted(found):
        try:
            sc = check_one(base, stem, found[stem])
            status[stem] = None
            by_id.setdefault(stem.rsplit("_", 1)[0], []).append((stem, sc))
        except Exception as e:      # one bad dataset must never hide the rest
            status[stem] = f"{type(e).__name__}: {e}"
    for entries in by_id.values():
        first_stem, first = entries[0]
        for stem, sc in entries[1:]:
            blocks = sidecar.invariant_diff(first, sc)
            if blocks:
                status[stem] = (f"disagrees with {first_stem} on the frequency-independent "
                                f"block(s) {', '.join(blocks)}")
    for stem in sorted(status):
        n_sets += 1
        if status[stem] is None:
            print(f"  PASS {base}/{stem}")
        else:
            n_failing += 1
            print(f"  FAIL {base}/{stem}: {status[stem]}")
print(f"{n_sets} dataset(s) checked, {n_failing} failing"
      + (f"; {n_bad_dirs} directory argument(s) failing" if n_bad_dirs else ""))
ok = n_failing == 0 and n_bad_dirs == 0
print("RESULT: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
