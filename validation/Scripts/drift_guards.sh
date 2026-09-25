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
NB_PY='import json, sys; nb = json.load(open(sys.argv[1])); sys.exit(1 if any(c.get("outputs") or c.get("execution_count") for c in nb["cells"] if c["cell_type"] == "code") else 0)'
report "notebooks without outputs" "$(git -C "$REPO" ls-files '*.ipynb' | while read -r nb; do
    $PY -c "$NB_PY" "$REPO/$nb" >/dev/null 2>&1 || echo "$nb"; done)"
exit $nfail
