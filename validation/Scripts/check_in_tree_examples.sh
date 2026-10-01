#!/usr/bin/env bash
# Configure and run examples against the source-tree library without an install.
# Pass any local CMake compiler/Geant4 options after <build-dir>. The dedicated
# scratch directory must be fresh or already use the Unix Makefiles generator.
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <repo> <build-dir> [cmake options...]" >&2
  exit 2
fi

repo=$1
build=$2
shift 2

cmake -G 'Unix Makefiles' -S "$repo" -B "$build" \
  -DBUILD_BBRSIM_EXAMPLES=ON \
  -DBUILD_BBRSIM_TOOLS=OFF \
  -DCMAKE_DISABLE_FIND_PACKAGE_BBRsim=TRUE \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON "$@"
cmake --build "$build" --target bbrsimTestWorld bbrsimLightPipe -j 4

python3 - "$repo" "$build/compile_commands.json" <<'PY'
import json
import pathlib
import sys

repo = pathlib.Path(sys.argv[1]).resolve()
commands = json.loads(pathlib.Path(sys.argv[2]).read_text())
for example in ('testworld', 'lightpipe'):
    entries = [entry for entry in commands
               if f'/examples/{example}/src/' in entry['file']]
    if not entries:
        raise SystemExit(f'no compile commands for {example}')
    for entry in entries:
        command = entry['command']
        if str(repo / 'library/include') not in command:
            raise SystemExit(f'{example} does not use source library headers')
        if str(repo / 'install/include/BBRsim') in command:
            raise SystemExit(f'{example} uses installed library headers')
PY

for example in testworld lightpipe; do
  link_command="$build/examples/$example/CMakeFiles/${example}Lib.dir/link.txt"
  if ! grep -qF '../../library/libBBRsim' "$link_command"; then
    echo "FAIL: $example does not link the source-tree BBRsim library" >&2
    exit 1
  fi
done

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
printf '/control/verbose 0\n' > "$scratch/smoke.mac"
env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH \
  "$build/examples/testworld/bbrsimTestWorld" "$scratch/smoke.mac" > "$scratch/testworld.log" 2>&1
env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH \
  "$build/examples/lightpipe/bbrsimLightPipe" "$scratch/smoke.mac" > "$scratch/lightpipe.log" 2>&1
env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH \
  bash "$repo/validation/Scripts/check_batch_cli.sh" \
  "$build/examples/testworld/bbrsimTestWorld" \
  "$build/examples/lightpipe/bbrsimLightPipe"

# Enabling source-tree examples must not change the top-level install surface.
cmake --build "$build" --target version > /dev/null
install_prefix=$(mktemp -d)
trap 'rm -rf "$scratch" "$install_prefix"' EXIT
if ! cmake --install "$build" --prefix "$install_prefix" > "$scratch/install.log" 2>&1; then
  cat "$scratch/install.log" >&2
  exit 1
fi
for artifact in bin/bbrsimTestWorld bin/bbrsimLightPipe \
                lib/libbbrsimTestWorld.dylib lib/libbbrsimLightPipe.dylib \
                lib/libbbrsimTestWorld.so lib/libbbrsimLightPipe.so \
                macros/bbrsimTestWorld macros/bbrsimLightPipe; do
  if [ -e "$install_prefix/$artifact" ]; then
    echo "FAIL: in-tree build installed example artifact $artifact" >&2
    exit 1
  fi
done
echo "PASS: in-tree examples build, link, and run without installed BBRsim"
