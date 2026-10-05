"""
check_dataset_sidecars.py
Every HFSS dataset (<id>_<freq>GHz_Ephi=0 directory) under the given waveguides
directories has a schema-1 sidecar <id>_<freq>GHz.dataset.json that passes
bbrsim.sidecar.check_full:
  * the checks BBRsim's C++ runs (F1-F11, F13; codes BBR024 / BBR025);
  * the full modes block re-derived from the cross-section;
  * the sha256, size and row count of the four CSVs;
  * the CSV checks C1-C5.
The frequencies of one ID must also agree on the frequency-independent blocks.

Usage:
    conda run -n bbrsim python validation/check_dataset_sidecars.py <waveguides-dir> [more ...]
"""
import os
import sys

from bbrsim import hfss, sidecar

dirs = sys.argv[1:] or [hfss.default_base_dir()]
n, fails, by_id = 0, [], {}
for base in dirs:
    stems = sorted(name[:-len("_Ephi=0")] for name in os.listdir(base)
                   if name.endswith("_Ephi=0") and os.path.isdir(os.path.join(base, name)))
    for stem in stems:
        n += 1
        dataset_id, label = stem.rsplit("_", 1)
        try:
            sc = sidecar.check_full(base, dataset_id, stem, hfss.parse_frequency_GHz(label) or -1.0)
            by_id.setdefault((base, dataset_id), []).append((stem, sc))
            print(f"  PASS {base}/{stem}")
        except (ValueError, OSError, KeyError) as e:
            fails.append(stem)
            print(f"  FAIL {base}/{stem}: {e}")
for (base, dataset_id), entries in sorted(by_id.items()):
    first_stem, first = entries[0]
    for stem, sc in entries[1:]:
        if not sidecar.same_invariant(first, sc):
            fails.append(stem)
            print(f"  FAIL {base}/{stem}: disagrees with {first_stem} on the frequency-independent physics")
print(f"{n} dataset(s) checked, {len(fails)} failing")
ok = n > 0 and not fails
print("RESULT: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
