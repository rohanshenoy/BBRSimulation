#!/bin/bash
# run_regression.sh — the one-command "did I break anything" for BBRsim.
#
#   validation/Scripts/run_regression.sh [BUILD_DIR]        (default: build)
#
# 1. Configures, builds and installs the library (BUILD_DIR, prefix BBR_PREFIX,
#    default <repo>/install), stops if the installed data or headers differ
#    from data/ or library/include (cmake --install never deletes), sources the
#    installed bbrsim_env.sh, then builds examples/testworld and
#    examples/lightpipe against it (BUILD_DIR/examples/<name>). Fails on build
#    errors; compiler warnings are a FAIL. The builds are incremental, so only
#    files compiled in this run are seen: use a fresh BUILD_DIR for a full
#    warnings audit. Release only: a BUILD_DIR already configured with another
#    build type (e.g. build-debug) is refused.
# 1b. Runs validation/Scripts/drift_guards.sh, checks that bbrsim imports from
#    the installed copy and that it matches tools/python/bbrsim, then builds
#    and runs validation/Scripts/consumer_smoke against the installed library.
# 2. Builds the mock HFSS tree when needed, checks that the real
#    data/waveguides holds only 500 GHz (leak guard), runs the eight cases in
#    parallel, each in BUILD_DIR/regression/<case>/ (BBRSIMDATA comes from the
#    installed env script): the five validation/G4Macros fixtures and three
#    example macros (reflectance.mac, planck.mac, lightpipe.mac; their command
#    lines are pinned by drift_guards.sh). Scans every log for GeomNav /
#    G4Exception / BBR0xx / LP002 messages, and checks that the frequency
#    case's BBR008 clamp warning fires once per side.
# 3. Runs every validation/check_*.py validator the fixtures feed and prints
#    one PASS/FAIL line per check.
#
# Exit code is the number of unexpected failures. A check listed in XFAIL is a
# known, documented red (see validation/README.md, PASS criteria). It is
# reported as XFAIL, without failing the run, only when its output shows the
# documented failure (see xfail_matches); any other failure, a crash included,
# is a FAIL, and an unexpected pass is flagged XPASS.
#
# Python is run as `conda run -n bbrsim python` (override with BBR_PYTHON).
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd -P)"
VAL="$REPO/validation"; VM="$VAL/G4Macros"
BUILD="${1:-build}"
case "$BUILD" in /*) ;; *) BUILD="$REPO/$BUILD" ;; esac
PY="${BBR_PYTHON:-conda run -n bbrsim python}"
REG="$BUILD/regression"
PREFIX="${BBR_PREFIX:-$REPO/install}"
case "$PREFIX" in /*) ;; *) PREFIX="$REPO/$PREFIX" ;; esac
# Canonical, as bbrsim_env.sh's BBRSIMINSTALL (pwd -P) and hence PYTHONPATH; this also
# removes symlinks (/tmp -> /private/tmp), .. and trailing slashes, which CMake collapses.
mkdir -p "$PREFIX" && PREFIX="$(cd "$PREFIX" >/dev/null && pwd -P)" || { echo "ERROR: cannot create $PREFIX"; exit 2; }
EXB="$BUILD/examples"
CLANG=(-DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release)
JOBS="${BBR_JOBS:-8}"

# Known red validators, "script|reason|signature" entries separated by ";", where
# signature is an ERE (no ";" or "|") that every failing row must match (see
# xfail_matches). BBR_XFAIL="" disables the list.
XFAIL="${BBR_XFAIL-check_cu_serov.py|HP_Cu alias RRR 6 is 13% low vs Serov under Drude (open decision, see validation/README.md)|^HP_Cu .*OUT OF TOLERANCE}"

fail=0; xfail=0; xpass=0; pass=0
line() { printf '%-6s %-30s %s\n' "$1" "$2" "$3"; }

if [ -f "$BUILD/CMakeCache.txt" ] && ! grep -q '^CMAKE_PROJECT_NAME:STATIC=BBRsim$' "$BUILD/CMakeCache.txt"; then
  echo "ERROR: $BUILD was configured by the pre-reorg single-project build; delete it and rerun."; exit 2
fi
# The configure below forces Release, which would silently turn a Debug build dir
# (clang-debug preset) into a Release one: refuse it instead.
if [ -f "$BUILD/CMakeCache.txt" ] && ! grep -qE '^CMAKE_BUILD_TYPE:[A-Z]*=Release$' "$BUILD/CMakeCache.txt"; then
  echo "ERROR: $BUILD is configured as CMAKE_BUILD_TYPE=$(sed -n 's/^CMAKE_BUILD_TYPE:[A-Z]*=//p' "$BUILD/CMakeCache.txt"); the runner builds Release only (default BUILD_DIR: build)."; exit 2
fi
echo "=== 1. build ($BUILD -> $PREFIX) ==="
clog="$(mktemp)"; blog="$(mktemp)"
step() {  # log label command... ; aborts the whole run on failure
  local log="$1" label="$2"; shift 2
  if ! "$@" >>"$log" 2>&1; then
    echo "BUILD FAILED ($label):"; grep -E "error|Error" "$log" | head -20; rm -f "$clog" "$blog"; exit 2
  fi
}
# BUILD_BBRSIM_TOOLS=ON: the validators import the installed bbrsim package, so a
# BUILD_DIR cached with it OFF must not skip reinstalling it.
step "$clog" "configure library" cmake -S "$REPO" -B "$BUILD" "${CLANG[@]}" -DCMAKE_INSTALL_PREFIX="$PREFIX" \
     -DBUILD_BBRSIM_TOOLS=ON
step "$blog" "build library"     cmake --build "$BUILD" -j"$JOBS"
step "$clog" "install library"   cmake --install "$BUILD"
# cmake --install never deletes, so a file once installed and since removed from
# the source survives in the prefix. The executables read the installed data copy
# (BBRSIMDATA from the env script): a stale dataset there (e.g. mock data written
# into data/ by mistake) would go unseen by the real-data leak guard below, which
# checks data/. The examples compile against the installed headers: one still
# including a removed or renamed header would build here and fail on a fresh
# install.
for pair in "data:share/BBRsim/data" "library/include:include/BBRsim"; do
  if ! diff -rq -x .DS_Store "$REPO/${pair%%:*}" "$PREFIX/${pair#*:}" >/dev/null 2>&1; then
    echo "ERROR: $PREFIX/${pair#*:} differs from $REPO/${pair%%:*} (cmake --install never deletes stale files); remove it and rerun."
    rm -f "$clog" "$blog"; exit 2
  fi
done
. "$PREFIX/share/BBRsim/bbrsim_env.sh" || { echo "ERROR: cannot source $PREFIX/share/BBRsim/bbrsim_env.sh"; exit 2; }
# Binaries resolve libBBRsim, Geant4 and their own example library through RPATH.
# dyld searches DYLD_LIBRARY_PATH by leaf name before @rpath, so the env script's
# $PREFIX/lib entry would let an example library previously installed there
# (cmake --install of an example; the presets use the same prefix) shadow the one
# just built in $EXB, and the fixtures would silently test stale example code.
# /bin/bash is SIP-protected, so this removes only what the env script added.
unset DYLD_LIBRARY_PATH LD_LIBRARY_PATH
# An example/consumer build dir configured against another prefix keeps its cached
# BBRsim_DIR (a new CMAKE_PREFIX_PATH is ignored), so it would keep the old prefix's
# headers, link line and RPATH: start it fresh.
for d in "$EXB/testworld" "$EXB/lightpipe" "$BUILD/consumer_smoke"; do
  if [ -f "$d/CMakeCache.txt" ] && ! grep -qxF "BBRsim_DIR:PATH=$PREFIX/lib/cmake/BBRsim" "$d/CMakeCache.txt"; then
    rm -rf "$d"
  fi
done
# NO_SYSTEM_FROM_IMPORTED passes imported include dirs as -I, not -isystem. As
# -isystem, BBRsim::BBRsim's would hide warnings in header-only library code that
# only the examples compile (BBRMaterials.hh, CADMesh.hh), so they would escape the
# zero-warnings gate. Geant4's dirs become -I too; its headers compile clean.
for ex in testworld lightpipe; do
  step "$clog" "configure $ex" cmake -S "$REPO/examples/$ex" -B "$EXB/$ex" "${CLANG[@]}" \
       -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_NO_SYSTEM_FROM_IMPORTED=ON
  step "$blog" "build $ex" cmake --build "$EXB/$ex" -j"$JOBS"
done
nwarn=$(grep -c -i "warning" "$blog" || true)
echo "build ok; compiler warnings: $nwarn"
[ "$nwarn" -eq 0 ] || { grep -i "warning" "$blog" | head -10; fail=$((fail+1)); line FAIL build "compiler warnings: $nwarn"; }
rm -f "$clog" "$blog"
TESTWORLD="$EXB/testworld/bbrsimTestWorld"
LIGHTPIPE="$EXB/lightpipe/bbrsimLightPipe"

echo "=== 1b. drift guards and consumer smoke ==="
gout=$("$VAL/Scripts/drift_guards.sh" 2>&1); grc=$?; echo "$gout"
gpass=$(echo "$gout" | grep -c '^PASS'); gfail=$(echo "$gout" | grep -c '^FAIL')
pass=$((pass + gpass)); fail=$((fail + gfail))
# Its exit code is its FAIL count. A script that cannot run (lost exec bit) or
# dies part-way would otherwise add no rows and leave the run green.
if [ $((gpass + gfail)) -eq 0 ] || [ "$grc" -ne "$gfail" ]; then
  line FAIL drift_guards "did not run cleanly (rc=$grc, $gpass PASS, $gfail FAIL)"; fail=$((fail+1))
fi
# With -c, Python puts the CWD first on the module search path, ahead of
# PYTHONPATH: a runner started in tools/python would import the source package,
# so run from /. And cmake --install never deletes, so a module removed or
# renamed in tools/python/bbrsim would survive in the prefix: compare the
# contents as well.
pyinst="$PREFIX/share/BBRsim/python/bbrsim"
if ! bp=$(cd / && $PY -c 'import bbrsim.paths, bbrsim.hfss, os; print(os.path.dirname(bbrsim.__file__))' 2>&1) \
   || [ "$bp" != "$pyinst" ]; then line FAIL "bbrsim importable" "got: $bp"; fail=$((fail+1))
elif ! diff -rq -x __pycache__ -x .DS_Store "$REPO/tools/python/bbrsim" "$pyinst" >/dev/null 2>&1; then
  line FAIL "bbrsim importable" "$pyinst differs from tools/python/bbrsim (cmake --install never deletes stale files); remove it and rerun"
  fail=$((fail+1))
else line PASS "bbrsim importable" "installed copy"; pass=$((pass+1)); fi
CS="$BUILD/consumer_smoke"; cslog="$BUILD/consumer_smoke.log"
if cmake -S "$VAL/Scripts/consumer_smoke" -B "$CS" "${CLANG[@]}" -DCMAKE_PREFIX_PATH="$PREFIX" >"$cslog" 2>&1 \
   && cmake --build "$CS" >>"$cslog" 2>&1 \
   && ( cd / && env -u BBRSIMDATA "$CS/consumer_smoke" ) 2>&1 | grep -qx "BBRsim data dir: $PREFIX/share/BBRsim/data"
then line PASS consumer_smoke "find_package(BBRsim) + link + compiled data default"; pass=$((pass+1))
else line FAIL consumer_smoke "see $cslog"; fail=$((fail+1)); fi

echo "=== 2. macros ==="
rm -rf "$REG"; mkdir -p "$REG"

# Mock multi-frequency HFSS tree for the frequency fixture. ~1.3 GB and ~5 s, so
# it is rebuilt only when it is missing or older than the real data or the
# generator. It lives in the build dir and is reached from a macro's run
# directory as ../mock_hfss via the symlink beside it.
MOCK="$BUILD/mock_hfss"
if [ ! -d "$MOCK/waveguides" ] || \
   [ -n "$(find "$REPO/data/waveguides" "$VAL/Scripts/make_mock_hfss_frequencies.py" \
            -newer "$MOCK/waveguides" -print -quit 2>/dev/null)" ]; then
  echo "generating mock HFSS frequency tree in $MOCK ..."
  # Written beside the final path and renamed only on success, so an interrupted
  # or failed generation cannot leave a partial tree that later runs would reuse.
  rm -rf "$MOCK" "$MOCK.tmp"
  if ! $PY "$VAL/Scripts/make_mock_hfss_frequencies.py" --src "$REPO/data/waveguides" \
         --dst "$MOCK.tmp" --ids InfParallelPlate_crack1Rohan InfParallelPlate_crack2 \
         >"$BUILD/mock_hfss.log" 2>&1 \
     || ! mv "$MOCK.tmp" "$MOCK"; then
    rm -rf "$MOCK.tmp"; echo "MOCK GENERATION FAILED:"; tail -5 "$BUILD/mock_hfss.log"; exit 2
  fi
fi
ln -s "$MOCK" "$REG/mock_hfss"

# The real tree must stay single-frequency: catches mock data written into data/.
# Run with -c: conda run does not forward stdin, so a heredoc would never execute
# and the guard could never fail. Run from / so the CWD cannot shadow the
# installed bbrsim (see the import check above).
LEAK_PY='import sys; from bbrsim import hfss; bad = [(i, g) for i in ("InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2") for g in [[f for f, _ in hfss.discover_frequencies(i, sys.argv[1])]] if g != [500.0]]; print(bad or "ok"); sys.exit(1 if bad else 0)'
if out=$(cd / && $PY -c "$LEAK_PY" "$REPO/data/waveguides" 2>&1); then
  line PASS "real-data leak guard" "data/waveguides holds only 500 GHz"; pass=$((pass+1))
else line FAIL "real-data leak guard" "$(echo "$out" | tail -1)"; fail=$((fail+1)); fi

run_macro() {  # case executable macro-path
  mkdir -p "$REG/$1" && ( cd "$REG/$1" && "$2" "$3" >run.log 2>&1; echo $? >exit.code )
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
# refl, planck and lp run the example macros themselves, pinned by the drift
# guard "regression macros pinned"; the other five are validation-only fixtures.
run_macro refl      "$TESTWORLD" "$REPO/examples/testworld/G4Macros/reflectance.mac" &
run_macro planck    "$TESTWORLD" "$REPO/examples/testworld/G4Macros/planck.mac"      &
run_macro wall      "$TESTWORLD" "$VM/Validation_CrackWall.mac"       &
run_macro exit      "$TESTWORLD" "$VM/Validation_WorldExit.mac"       &
run_macro transmit  "$TESTWORLD" "$VM/Validation_CrackTransmit.mac"   &
run_macro oblique   "$TESTWORLD" "$VM/Validation_CrackOblique.mac"    &
run_macro frequency "$TESTWORLD" "$VM/Validation_CrackFrequency.mac"  &
run_macro lp        "$LIGHTPIPE" "$REPO/examples/lightpipe/G4Macros/lightpipe.mac"   &
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
# of the frequency fixture). crack2 may add at most one of each; the tokens carry
# the dataset and side so the count cannot depend on which run clamps first.
for side in low high; do
  n=$(grep -c "BBR008 dataset=InfParallelPlate_crack1Rohan side=$side" "$REG/frequency/run.log" || true)
  if [ "$n" -eq 1 ]; then line PASS "clamp warning $side" "[frequency] BBR008 once"; pass=$((pass+1))
  else line FAIL "clamp warning $side" "[frequency] BBR008 count $n, expected 1"; fail=$((fail+1)); fi
done

echo "=== 3. validators ==="
# An XFAIL holds only for the documented failure: the validator reached its
# "RESULT: FAIL" line, and every failing row (a FAIL or OUT OF TOLERANCE line
# other than RESULT) matches the entry's signature, at least one of them. A
# crash, an import error or another row going out of tolerance is a plain FAIL.
xfail_matches() {  # output signature-ERE (empty: the RESULT line alone)
  local rows
  printf '%s\n' "$1" | grep -qE '^ *RESULT *: *FAIL' || return 1
  [ -n "$2" ] || return 0
  rows=$(printf '%s\n' "$1" | grep -E 'FAIL|OUT OF TOLERANCE' | grep -vE '^ *RESULT')
  [ -n "$rows" ] && ! printf '%s\n' "$rows" | grep -qvE "$2"
}
check() {  # case[/file] script [extra args...]; "-" = no ROOT input
  # The input is the case's output/bbr.root, or output/<file> when one is named.
  local d="$1" script="$2"; shift 2
  local out res f="bbr.root"
  case "$d" in */*) f="${d#*/}" ;; esac
  if [ "$d" = "-" ]; then out=$($PY "$VAL/$script" "$@" 2>&1)
  else                    out=$($PY "$VAL/$script" "$REG/${d%%/*}/output/$f" "$@" 2>&1); fi
  local rc=$?
  res=$(echo "$out" | grep -E "RESULT|PASS|FAIL" | tail -1 | sed 's/^ *//')
  local entry reason sig
  entry=$(printf '%s\n' "$XFAIL" | tr ';' '\n' | grep "^$script|" | head -1)
  reason=$(echo "$entry" | cut -d'|' -f2); sig=$(echo "$entry" | cut -d'|' -f3-)
  if [ $rc -eq 0 ]; then
    if [ -n "$entry" ]; then line XPASS "$script" "[$d] expected to fail but passed — remove it from XFAIL"; xpass=$((xpass+1)); fail=$((fail+1))
    else line PASS "$script" "[$d] $res"; pass=$((pass+1)); fi
  elif [ -n "$entry" ] && xfail_matches "$out" "$sig"; then
    line XFAIL "$script" "[$d] $reason"; xfail=$((xfail+1))
  else
    line FAIL "$script" "[$d] $res${entry:+ (not the documented XFAIL)}"; fail=$((fail+1))
    echo "$out" | tail -6 | sed 's/^/       /'
  fi
}
cd "$REG" || exit 2   # validators write plots into the CWD; keep them out of the repo
# check_reflectance takes --root rather than a positional path, so it is run directly.
out=$($PY "$VAL/check_reflectance.py" --root "$REG/refl/output/bbr.root" 2>&1); rc=$?
res=$(echo "$out" | grep -E "PASS|FAIL" | tail -1 | sed 's/^ *//')
if [ $rc -eq 0 ]; then line PASS check_reflectance.py "[refl] $res"; pass=$((pass+1)); else line FAIL check_reflectance.py "[refl] $res"; fail=$((fail+1)); fi
check refl     check_invariants.py
check planck   check_planck_spectrum.py --temp 4
check planck   check_nreflect.py
check planck   check_angle_distribution.py
check planck   check_invariants.py
check wall     check_crack_wall_reflection.py
check wall     check_invariants.py
# the world-exit fixture crosses no boundary by design, so its file holds no
# crossings; the flag waives only that requirement of the metal invariant.
check exit     check_invariants.py --allow-no-crossings
check transmit check_crack_transmittance.py
check transmit check_invariants.py
# the transmit fixture's second run (Planck, both cracks) has its own file.
check transmit/bbr_ratio.root check_crack_ratio.py
check lp       check_invariants.py
# the oblique fixture writes one file per run (output/bbr_oblique_rNN.root);
# check_crack_oblique reads the whole directory, check_invariants runs on every
# per-run file.
out=$($PY "$VAL/check_crack_oblique.py" "$REG/oblique/output" 2>&1); rc=$?
res=$(echo "$out" | grep -E "^RESULT" | tail -1)
if [ $rc -eq 0 ]; then line PASS check_crack_oblique.py "[oblique] $res"; pass=$((pass+1))
else line FAIL check_crack_oblique.py "[oblique] $res"; fail=$((fail+1)); echo "$out" | grep -E "^  FAIL" | head -8 | sed 's/^/       /'; fi
obl_bad=0; obl_n=0
for f in "$REG"/oblique/output/bbr_oblique_r*.root; do
  [ -f "$f" ] || continue
  obl_n=$((obl_n+1))
  out=$($PY "$VAL/check_invariants.py" "$f" 2>&1) || { obl_bad=$((obl_bad+1)); line FAIL check_invariants.py "[oblique/$(basename "$f")] $(echo "$out" | grep -E "^RESULT" | tail -1)"; }
done
if [ "$obl_n" -gt 0 ] && [ "$obl_bad" -eq 0 ]; then line PASS "oblique invariants" "[oblique] check_invariants on $obl_n per-run files"; pass=$((pass+1))
else fail=$((fail+obl_bad)); [ "$obl_n" -eq 0 ] && { line FAIL "oblique invariants" "[oblique] no per-run files found"; fail=$((fail+1)); }; fi
out=$($PY "$VAL/check_crack_frequency.py" "$REG/frequency/output" --data-dir "$REG/mock_hfss" 2>&1); rc=$?
res=$(echo "$out" | grep -E "^RESULT" | tail -1)
if [ $rc -eq 0 ]; then line PASS check_crack_frequency.py "[frequency] $res"; pass=$((pass+1))
else line FAIL check_crack_frequency.py "[frequency] $res"; fail=$((fail+1)); echo "$out" | grep -E "^  FAIL" | head -8 | sed 's/^/       /'; fi
frq_bad=0; frq_n=0
for f in "$REG"/frequency/output/bbr_freq_r*.root; do
  [ -f "$f" ] || continue
  frq_n=$((frq_n+1))
  out=$($PY "$VAL/check_invariants.py" "$f" 2>&1) || { frq_bad=$((frq_bad+1)); line FAIL check_invariants.py "[frequency/$(basename "$f")] $(echo "$out" | grep -E "^RESULT" | tail -1)"; }
done
if [ "$frq_n" -gt 0 ] && [ "$frq_bad" -eq 0 ]; then line PASS "frequency invariants" "[frequency] check_invariants on $frq_n per-run files"; pass=$((pass+1))
else fail=$((fail+frq_bad)); [ "$frq_n" -eq 0 ] && { line FAIL "frequency invariants" "[frequency] no per-run files found"; fail=$((fail+1)); }; fi
check -        check_physics.py
check -        check_cu_serov.py

echo "=== summary ==="
echo "pass=$pass  fail=$fail  xfail=$xfail  xpass=$xpass   (outputs under $REG)"
exit $fail
