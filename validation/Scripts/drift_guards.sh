#!/bin/bash
# drift_guards.sh — source-tree consistency checks that no compiler or validator
# catches. Prints one "PASS|FAIL <label> <detail>" line per guard; the exit code
# is the number of FAIL lines. Called by run_regression.sh; runnable alone.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd -P)"
PY="${BBR_PYTHON:-conda run -n bbrsim python}"
nfail=0
report() {  # label missing-list
  if [ -z "$2" ]; then printf '%-6s %-30s %s\n' PASS "$1" ok
  else printf '%-6s %-30s %s\n' FAIL "$1" "$(echo "$2" | head -3 | tr '\n' ' ')"; nfail=$((nfail+1)); fi
}
unlisted() {  # cmakelists files... : files whose name the CMakeLists does not mention
  local cm="$1"; shift
  for f in "$@"; do
    [ -e "$f" ] || continue
    b="$(basename "$f")"
    grep -qF "/$b" "$cm" || grep -qwF "$b" "$cm" || echo "${f#$REPO/}"
  done
}
report "library file lists" "$(unlisted "$REPO/library/CMakeLists.txt" "$REPO"/library/include/* "$REPO"/library/src/*)"
report "example file lists" "$(for ex in testworld lightpipe; do
    unlisted "$REPO/examples/$ex/CMakeLists.txt" "$REPO/examples/$ex"/src/*.cc "$REPO/examples/$ex"/G4Macros/*.mac; done)"
report "library names no example" "$(grep -lE '(TestWorld|LightPipe)(DetectorConstruction|ActionInitialization|PrimaryGeneratorAction|SteppingAction|RunAction|Messenger)|BBR(Test|RunAction|LightPipe)' \
    "$REPO"/library/include/* "$REPO"/library/src/* 2>/dev/null | sed "s#^$REPO/##")"
report "every fixture is run" "$(for f in "$REPO"/validation/G4Macros/Validation_*.mac; do
    grep -qF "$(basename "$f")" "$REPO/validation/Scripts/run_regression.sh" || basename "$f"; done)"
# The runner runs these four example macros directly (cases refl, planck,
# config_mt, lp), and their validators rely on the seeds, event counts and
# settings. Each pin is the sha256 of the macro's lines with full-line comments
# and blank lines stripped, so only edits to those lines pass; any other change
# fails, including an inline comment or trailing whitespace on a command line.
# The guard also fails if the runner no longer runs the macro. An installed copy
# of validation/ has no examples/ unless INSTALL_EXAMPLES was on, so there it is
# skipped.
PINS="examples/testworld/G4Macros/reflectance.mac cf4e27f513d8eccd9dbf28f0834f88c1d83b06f497de0a35169915bb2d736317
examples/testworld/G4Macros/planck.mac      ba261b402f09768d9fbdc464b177ec285354e2f9d733e64da91fecad5623b6d5
examples/testworld/G4Macros/config_mt.mac   36dda6c6b9699527adc6e84e9880ad01fb662e4f20492a254ef5636876baf57c
examples/lightpipe/G4Macros/lightpipe.mac   e2fbb3c3df2ba03a604c96356e2a69e9f675327dc0b664c1c8120c1d570b3b68"
sha256() { if command -v sha256sum >/dev/null 2>&1; then sha256sum; else shasum -a 256; fi | cut -d' ' -f1; }
if [ -d "$REPO/examples" ]; then
  report "regression macros pinned" "$(printf '%s\n' "$PINS" | while read -r m pin; do
    grep -F "\"\$REPO/$m\"" "$REPO/validation/Scripts/run_regression.sh" | grep -q '^run_macro ' ||
      echo "$m is not run by run_regression.sh: it is a regression input;"
    if [ ! -f "$REPO/$m" ]; then echo "$m is missing: it is a regression input (run_regression.sh);"; continue; fi
    [ "$(grep -vE '^[[:space:]]*(#|$)' "$REPO/$m" | sha256)" = "$pin" ] ||
      echo "$m changed: it is a regression input; change it only together with its validators, then update its pin in validation/Scripts/drift_guards.sh;"
    done)"
else
  printf '%-6s %-30s %s\n' PASS "regression macros pinned" "skipped: no examples/ here (installed copy)"
fi
NB_PY='import json, sys; nb = json.load(open(sys.argv[1])); sys.exit(1 if any(c.get("outputs") or c.get("execution_count") for c in nb["cells"] if c["cell_type"] == "code") else 0)'
report "notebooks without outputs" "$(git -C "$REPO" ls-files '*.ipynb' | while read -r nb; do
    $PY -c "$NB_PY" "$REPO/$nb" >/dev/null 2>&1 || echo "$nb"; done)"
report "tools file lists" "$(unlisted "$REPO/tools/CMakeLists.txt" "$REPO"/tools/python/bbrsim/*.py "$REPO"/tools/plot_*.py)"
exit $nfail
