#!/bin/bash
# run_regression.sh — the one-command "did I break anything" for BBRsim.
#
#   scripts/run_regression.sh [BUILD_DIR]        (default: build)
#
# 1. Incremental build of BBRSim + BBRLightPipe in BUILD_DIR (fails on errors,
#    reports compiler warnings).
# 2. Runs the regression macros in parallel, each in its own directory under
#    BUILD_DIR/regression/<name>/ (BBRSIMDATA defaults to <repo>/data, so the
#    executables and the Python validators read the same tree), and scans every
#    log for GeomNav / G4Exception / BBR00x / LP002 messages.
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
# The executables' compiled-in data default is the INSTALLED <prefix>/share/BBRsim/data,
# and this runner builds without installing. Unless BBRSIMDATA is already set, point
# the executables and the Python validators at this checkout's data/.
export BBRSIMDATA="${BBRSIMDATA:-$REPO/data}"

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

# Mock multi-frequency HFSS tree for crack_frequency.mac. ~1.3 GB and ~5 s, so
# it is rebuilt only when it is missing or older than the real data or the
# generator. It lives in the build dir and is reached from a macro's run
# directory as ../mock_hfss via the symlink beside it.
MOCK="$BUILD/mock_hfss"
if [ ! -d "$MOCK/waveguides" ] || \
   [ -n "$(find "$REPO/data/waveguides" "$REPO/scripts/make_mock_hfss_frequencies.py" \
            -newer "$MOCK/waveguides" -print -quit 2>/dev/null)" ]; then
  echo "generating mock HFSS frequency tree in $MOCK ..."
  rm -rf "$MOCK"
  $PY "$REPO/scripts/make_mock_hfss_frequencies.py" --src "$REPO/data/waveguides" \
      --dst "$MOCK" --ids InfParallelPlate_crack1Rohan InfParallelPlate_crack2 \
      >"$BUILD/mock_hfss.log" 2>&1 \
    || { echo "MOCK GENERATION FAILED:"; tail -5 "$BUILD/mock_hfss.log"; exit 2; }
fi
ln -s "$MOCK" "$REG/mock_hfss"

# The real tree must stay single-frequency: catches mock data written into data/.
# Run with -c: conda run does not forward stdin, so a heredoc would never execute
# and the guard could never fail.
LEAK_PY='import sys; sys.path.insert(0, sys.argv[1] + "/analysis"); from bbrsim import hfss; bad = [(i, g) for i in ("InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2") for g in [[f for f, _ in hfss.discover_frequencies(i, sys.argv[2])]] if g != [500.0]]; print(bad or "ok"); sys.exit(1 if bad else 0)'
if out=$($PY -c "$LEAK_PY" "$REPO" "$REPO/data/waveguides" 2>&1); then
  line PASS "real-data leak guard" "data/waveguides holds only 500 GHz"; pass=$((pass+1))
else line FAIL "real-data leak guard" "$(echo "$out" | tail -1)"; fail=$((fail+1)); fi

run_macro() {  # name executable macro
  mkdir -p "$REG/$1" && ( cd "$REG/$1" && "$BUILD/$2" "$REPO/$3" >run.log 2>&1; echo $? >exit.code )
}

# Count log lines that indicate trouble. A tolerated exception code (the
# expected BBR008 clamp warnings of the frequency case) is subtracted from both
# the code hits and the two-line G4Exception banner, which a single grep cannot
# do: a JustWarning prints "G4Exception-START" and "*** G4Exception : CODE" on
# separate lines, and threads interleave.
scan_log() {  # case-dir [tolerated-code]
  local log="$REG/$1/run.log" tol="${2:-}"
  local n_geom n_lp n_start n_bbr n_tol=0
  n_geom=$(grep -c -i "GeomNav" "$log" || true)
  n_lp=$(grep -c "LP002" "$log" || true)
  n_start=$(grep -c "G4Exception-START" "$log" || true)
  if [ -n "$tol" ]; then
    n_bbr=$(grep -E "BBR0[0-9][0-9]" "$log" | grep -v -c "$tol" || true)
    n_tol=$(grep -c -E "G4Exception : $tol" "$log" || true)
  else
    n_bbr=$(grep -c -E "BBR0[0-9][0-9]" "$log" || true)
  fi
  echo $(( n_geom + n_lp + (n_start - n_tol) + n_bbr ))
}
run_macro refl     BBRSim       reflectance.mac    &
run_macro planck   BBRSim       planck.mac         &
run_macro wall     BBRSim       crack_wall.mac     &
run_macro exit     BBRSim       world_exit.mac     &
run_macro transmit BBRSim       crack_transmit.mac &
run_macro oblique  BBRSim       crack_oblique.mac  &
run_macro frequency BBRSim      crack_frequency.mac &
run_macro lp       BBRLightPipe lightpipe.mac      &
wait
for d in refl planck wall exit transmit oblique frequency lp; do
  code=$(cat "$REG/$d/exit.code")
  tol=""; [ "$d" = "frequency" ] && tol="BBR008"   # expected clamp warnings
  nbad=$(scan_log "$d" "$tol")
  if [ "$code" -eq 0 ] && [ "$nbad" -eq 0 ] && ls "$REG/$d"/output/*.root >/dev/null 2>&1; then
    line PASS "run:$d" "exit 0, no GeomNav/G4Exception/BBR0xx/LP002${tol:+ (except $tol)}"; pass=$((pass+1))
  else
    line FAIL "run:$d" "exit $code, flagged log lines: $nbad"; fail=$((fail+1))
    grep -E "GeomNav|G4Exception : |BBR0[0-9][0-9]|LP002" "$REG/$d/run.log" | head -3
  fi
done
# The clamp warning must fire exactly once per side for crack1 (runs 13 and 14
# of crack_frequency.mac). crack2 may add at most one of each; the tokens carry
# the dataset and side so the count cannot depend on which run clamps first.
for side in low high; do
  n=$(grep -c "BBR008 dataset=InfParallelPlate_crack1Rohan side=$side" "$REG/frequency/run.log" || true)
  if [ "$n" -eq 1 ]; then line PASS "clamp warning $side" "[frequency] BBR008 once"; pass=$((pass+1))
  else line FAIL "clamp warning $side" "[frequency] BBR008 count $n, expected 1"; fail=$((fail+1)); fi
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
# crack_oblique.mac writes one file per run (output/bbr_oblique_rNN.root);
# check_crack_oblique reads the whole directory, the two invariant checks run
# on every per-run file.
out=$($PY "$REPO/scripts/check_crack_oblique.py" "$REG/oblique/output" 2>&1); rc=$?
res=$(echo "$out" | grep -E "^RESULT" | tail -1)
if [ $rc -eq 0 ]; then line PASS check_crack_oblique.py "[oblique] $res"; pass=$((pass+1))
else line FAIL check_crack_oblique.py "[oblique] $res"; fail=$((fail+1)); echo "$out" | grep -E "^  FAIL" | head -8 | sed 's/^/       /'; fi
obl_bad=0; obl_n=0
for f in "$REG"/oblique/output/bbr_oblique_r*.root; do
  [ -f "$f" ] || continue
  obl_n=$((obl_n+1))
  for s in check_no_photons_in_metal.py check_term_status.py; do
    out=$($PY "$REPO/scripts/$s" "$f" 2>&1) || { obl_bad=$((obl_bad+1)); line FAIL "$s" "[oblique/$(basename "$f")] $(echo "$out" | grep -E "^RESULT" | tail -1)"; }
  done
done
if [ "$obl_n" -gt 0 ] && [ "$obl_bad" -eq 0 ]; then line PASS "oblique invariants" "[oblique] no_photons_in_metal + term_status on $obl_n per-run files"; pass=$((pass+1))
else fail=$((fail+obl_bad)); [ "$obl_n" -eq 0 ] && { line FAIL "oblique invariants" "[oblique] no per-run files found"; fail=$((fail+1)); }; fi
out=$($PY "$REPO/scripts/check_crack_frequency.py" "$REG/frequency/output" --data-dir "$REG/mock_hfss" 2>&1); rc=$?
res=$(echo "$out" | grep -E "^RESULT" | tail -1)
if [ $rc -eq 0 ]; then line PASS check_crack_frequency.py "[frequency] $res"; pass=$((pass+1))
else line FAIL check_crack_frequency.py "[frequency] $res"; fail=$((fail+1)); echo "$out" | grep -E "^  FAIL" | head -8 | sed 's/^/       /'; fi
frq_bad=0; frq_n=0
for f in "$REG"/frequency/output/bbr_freq_r*.root; do
  [ -f "$f" ] || continue
  frq_n=$((frq_n+1))
  for s in check_no_photons_in_metal.py check_term_status.py; do
    out=$($PY "$REPO/scripts/$s" "$f" 2>&1) || { frq_bad=$((frq_bad+1)); line FAIL "$s" "[frequency/$(basename "$f")] $(echo "$out" | grep -E "^RESULT" | tail -1)"; }
  done
done
if [ "$frq_n" -gt 0 ] && [ "$frq_bad" -eq 0 ]; then line PASS "frequency invariants" "[frequency] no_photons_in_metal + term_status on $frq_n per-run files"; pass=$((pass+1))
else fail=$((fail+frq_bad)); [ "$frq_n" -eq 0 ] && { line FAIL "frequency invariants" "[frequency] no per-run files found"; fail=$((fail+1)); }; fi
check -        check_physics.py
check -        check_cu_serov.py

echo "=== summary ==="
echo "pass=$pass  fail=$fail  xfail=$xfail  xpass=$xpass   (outputs under $REG)"
exit $fail
