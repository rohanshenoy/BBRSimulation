#!/bin/bash
# run_regression.sh — the one-command "did I break anything" for BBRsim.
#
#   scripts/run_regression.sh [BUILD_DIR]        (default: build)
#
# 1. Incremental build of BBRSim + BBRLightPipe in BUILD_DIR (fails on errors,
#    reports compiler warnings).
# 2. Runs the regression macros in parallel, each in its own directory under
#    BUILD_DIR/regression/<name>/ (a `data` symlink beside them satisfies the
#    executables' ../data/waveguides lookup), and scans every log for
#    GeomNav / G4Exception / BBR00x / LP002 messages.
# 3. Runs every scripts/check_*.py validator against the output it belongs to
#    and prints one PASS/FAIL line per check.
#
# Exit code is the number of unexpected failures. A check listed in XFAIL is a
# known, documented red (the open HP_Cu decision); it is reported but
# does not fail the run, and is flagged XPASS if it unexpectedly passes.
#
# Python is run as `conda run -n bbrsim python` (override with BBR_PYTHON).
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd -P)"
BUILD="${1:-build}"
case "$BUILD" in /*) ;; *) BUILD="$REPO/$BUILD" ;; esac
PY="${BBR_PYTHON:-conda run -n bbrsim python}"
REG="$BUILD/regression"
JOBS="${BBR_JOBS:-8}"

# Known red validators, "script|reason" entries separated by ";". BBR_XFAIL="" disables the list.
XFAIL="${BBR_XFAIL-check_cu_serov.py|HP_Cu alias RRR 6 is 13% low vs Serov under Drude (open decision)}"

fail=0; xfail=0; xpass=0; pass=0
line() { printf '%-6s %-30s %s\n' "$1" "$2" "$3"; }

[ -f "$BUILD/CMakeCache.txt" ] || { echo "ERROR: $BUILD is not a configured build directory (run cmake --preset clang-release first)"; exit 2; }

echo "=== 1. build ($BUILD) ==="
blog="$(mktemp)"
if ! make -C "$BUILD" -j"$JOBS" >"$blog" 2>&1; then
  echo "BUILD FAILED:"; grep -E "error|Error" "$blog" | head -20; rm -f "$blog"; exit 2
fi
nwarn=$(grep -c -i "warning" "$blog" || true)
echo "build ok; compiler warnings: $nwarn"
[ "$nwarn" -eq 0 ] || { grep -i "warning" "$blog" | head -10; fail=$((fail+1)); line FAIL build "compiler warnings: $nwarn"; }
rm -f "$blog"

echo "=== 2. macros ==="
rm -rf "$REG"; mkdir -p "$REG"; ln -s "$REPO/data" "$REG/data"
run_macro() {  # name executable macro
  mkdir -p "$REG/$1" && ( cd "$REG/$1" && "$BUILD/$2" "$REPO/$3" >run.log 2>&1; echo $? >exit.code )
}
run_macro refl     BBRSim       reflectance.mac    &
run_macro planck   BBRSim       planck.mac         &
run_macro wall     BBRSim       crack_wall.mac     &
run_macro exit     BBRSim       world_exit.mac     &
run_macro transmit BBRSim       crack_transmit.mac &
run_macro lp       BBRLightPipe lightpipe.mac      &
wait
for d in refl planck wall exit transmit lp; do
  code=$(cat "$REG/$d/exit.code"); nbad=$(grep -c -i "GeomNav\|G4Exception-START\|BBR00\|LP002" "$REG/$d/run.log" || true)
  if [ "$code" -eq 0 ] && [ "$nbad" -eq 0 ] && [ -f "$REG/$d/output/bbr.root" ]; then
    line PASS "run:$d" "exit 0, no GeomNav/G4Exception/BBR00x/LP002"; pass=$((pass+1))
  else
    line FAIL "run:$d" "exit $code, flagged log lines: $nbad"; fail=$((fail+1))
    grep -i "GeomNav\|G4Exception-START\|BBR00\|LP002" "$REG/$d/run.log" | head -3
  fi
done

echo "=== 3. validators ==="
check() {  # macro-dir script [extra args...]; macro-dir "-" = no ROOT input
  local d="$1" script="$2"; shift 2
  local out res
  if [ "$d" = "-" ]; then out=$($PY "$REPO/scripts/$script" "$@" 2>&1)
  else                    out=$($PY "$REPO/scripts/$script" "$REG/$d/output/bbr.root" "$@" 2>&1); fi
  local rc=$?
  res=$(echo "$out" | grep -E "RESULT|PASS|FAIL" | tail -1 | sed 's/^ *//')
  local reason; reason=$(printf '%s\n' "$XFAIL" | tr ';' '\n' | grep "^$script|" | cut -d'|' -f2-)
  if [ $rc -eq 0 ]; then
    if [ -n "$reason" ]; then line XPASS "$script" "[$d] expected to fail but passed — remove it from XFAIL"; xpass=$((xpass+1)); fail=$((fail+1))
    else line PASS "$script" "[$d] $res"; pass=$((pass+1)); fi
  else
    if [ -n "$reason" ]; then line XFAIL "$script" "[$d] $reason"; xfail=$((xfail+1))
    else line FAIL "$script" "[$d] $res"; fail=$((fail+1)); echo "$out" | tail -6 | sed 's/^/       /'; fi
  fi
}
cd "$REG" || exit 2   # validators write plots into the CWD; keep them out of the repo
# check_reflectance takes --root rather than a positional path, so it is run directly.
out=$($PY "$REPO/scripts/check_reflectance.py" --root "$REG/refl/output/bbr.root" 2>&1); rc=$?
res=$(echo "$out" | grep -E "PASS|FAIL" | tail -1 | sed 's/^ *//')
if [ $rc -eq 0 ]; then line PASS check_reflectance.py "[refl] $res"; pass=$((pass+1)); else line FAIL check_reflectance.py "[refl] $res"; fail=$((fail+1)); fi
check refl     check_term_status.py
check refl     check_no_photons_in_metal.py
check planck   check_planck_spectrum.py --temp 4
check planck   check_nreflect.py
check planck   check_angle_distribution.py
check planck   check_term_status.py
check planck   check_no_photons_in_metal.py
check wall     check_crack_wall_reflection.py
check wall     check_no_photons_in_metal.py
check wall     check_term_status.py
check exit     check_term_status.py
check transmit check_crack_transmittance.py
check transmit check_no_photons_in_metal.py
check transmit check_term_status.py
check lp       check_no_photons_in_metal.py
check lp       check_term_status.py
check -        check_physics.py
check -        check_cu_serov.py

echo "=== summary ==="
echo "pass=$pass  fail=$fail  xfail=$xfail  xpass=$xpass   (outputs under $REG)"
exit $fail
