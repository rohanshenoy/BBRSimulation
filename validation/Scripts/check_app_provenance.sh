#!/usr/bin/env bash
# A standalone example must refresh its own source identity after a local edit.
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "usage: $0 <repo> <installed-prefix> <scratch-dir> [cmake options...]" >&2
  exit 2
fi

repo=$1
prefix=$2
scratch=$3
shift 3
copy="$scratch/example"
build="$scratch/build"
library=""
for candidate in "$prefix"/lib/libBBRsim.*; do
  if [ -f "$candidate" ]; then library=$candidate; break; fi
done
if [ -z "$library" ]; then
  echo "FAIL: installed BBRsim library missing from $prefix" >&2
  exit 1
fi
library_hash_before=$(python3 - "$library" <<'PY'
import hashlib
import pathlib
import sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
mkdir -p "$copy"
cp "$repo/examples/testworld/CMakeLists.txt" "$copy/"
cp "$repo/examples/testworld/bbrsimTestWorld.cc" "$copy/"
cp -R "$repo/examples/testworld/include" "$repo/examples/testworld/src" \
      "$repo/examples/testworld/G4Macros" "$copy/"

env -u BBRSIMINSTALL cmake -G 'Unix Makefiles' -S "$copy" -B "$build" \
  -DCMAKE_PREFIX_PATH="$prefix" -DCMAKE_INSTALL_PREFIX="$prefix" "$@"
cmake --build "$build" -j 4

header="$build/BBRApplicationInfo.hh"
if [ ! -f "$header" ]; then
  echo "FAIL: standalone example has no generated application provenance" >&2
  exit 1
fi
first=$(rg '^#define BBRSIM_APP_SOURCE_FINGERPRINT ' "$header")
if [ -z "$first" ]; then
  echo "FAIL: application fingerprint is empty" >&2
  exit 1
fi
printf '/run/numberOfThreads 1\n/run/initialize\n/bbr/gun/mode true\n/bbr/gun/posY 1\n/bbr/gun/posZ 1\n/run/beamOn 1\n' > "$scratch/run.mac"
mkdir -p "$scratch/before" "$scratch/after"
(cd "$scratch/before" && env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH \
  BBRSIMDATA="$prefix/share/BBRsim/data" G4FORCENUMBEROFTHREADS=1 \
  "$build/bbrsimTestWorld" "$scratch/run.mac" > "$scratch/before.log" 2>&1)

printf '\n// Provenance fixture: local source changed.\n' >> "$copy/bbrsimTestWorld.cc"
cmake --build "$build" -j 4
second=$(rg '^#define BBRSIM_APP_SOURCE_FINGERPRINT ' "$header")
if [ "$first" = "$second" ]; then
  echo "FAIL: standalone source edit did not refresh application fingerprint" >&2
  exit 1
fi
(cd "$scratch/after" && env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH \
  BBRSIMDATA="$prefix/share/BBRsim/data" G4FORCENUMBEROFTHREADS=1 \
  "$build/bbrsimTestWorld" "$scratch/run.mac" > "$scratch/after.log" 2>&1)
python3 - "$scratch/before/output/bbr.metadata.json" \
          "$scratch/after/output/bbr.metadata.json" "$header" <<'PY'
import json
import pathlib
import re
import sys

before = json.loads(pathlib.Path(sys.argv[1]).read_text())['build']
after = json.loads(pathlib.Path(sys.argv[2]).read_text())['build']
header = pathlib.Path(sys.argv[3]).read_text()
match = re.search(r'BBRSIM_APP_SOURCE_FINGERPRINT "([^"]+)"', header)
if not match or before['application_source_fingerprint'] == after['application_source_fingerprint']:
    raise SystemExit('application metadata did not change after the source edit')
if after['application_source_fingerprint'] != match.group(1):
    raise SystemExit('application metadata does not match the rebuilt header')
if before['source_fingerprint'] != after['source_fingerprint']:
    raise SystemExit('library provenance changed during an example-only rebuild')
PY
library_hash_after=$(python3 - "$library" <<'PY'
import hashlib
import pathlib
import sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
if [ "$library_hash_before" != "$library_hash_after" ]; then
  echo "FAIL: standalone rebuild changed the installed BBRsim library" >&2
  exit 1
fi

echo "PASS: standalone application provenance follows source edits"
