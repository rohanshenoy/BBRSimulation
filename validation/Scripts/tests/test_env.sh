#!/bin/bash
# test_env.sh REPO PREFIX — env-script self-location in bash, zsh, tcsh and dash:
# 22 cases against the source tree REPO and the install PREFIX. Prints one
# "ok|FAIL|SKIP <case>" line per case (SKIP: that case's shell is not
# installed), then "ENV OK" or "ENV BAD"; the exit code is the FAIL count.
# Scratch files go to a temporary directory under ${TMPDIR:-/tmp}.
R="$1"; P="$2"; I="$P/share/BBRsim"; bad=0
SRC="$R|$R/build/library|$R/library/include|$R/data|$R/tools/python"
INS="$I|$P/lib|$P/include/BBRsim|$I/data|$I/python"
Q='echo "$BBRSIMINSTALL|$BBRSIMLIB|$BBRSIMINCLUDE|$BBRSIMDATA|$PYTHONPATH"'
try_in() {  # dir label expected command...  (the first bash/zsh/tcsh/dash word is the shell)
  local got a
  for a in "${@:4}"; do
    case "$a" in
      bash|zsh|tcsh|dash)
        command -v "$a" >/dev/null 2>&1 || { echo "SKIP $2 ($a not installed)"; return; }
        break ;;
    esac
  done
  got=$(cd "$1" && env -u BBRSIMINSTALL -u PYTHONPATH "${@:4}" 2>&1 | tail -1)
  if [ "$got" = "$3" ]; then echo "ok   $2"; else echo "FAIL $2: got '$got'"; bad=$((bad+1)); fi
}
try() { try_in / "$@"; }  # label expected command...
try "bash source-tree"  "$SRC" bash -c ". '$R/bbrsim_env.sh' && $Q"
try "bash install-tree" "$INS" bash -c ". '$I/bbrsim_env.sh' && $Q"
try "zsh source-tree"   "$SRC" zsh -f -c ". '$R/bbrsim_env.sh' && $Q"
try "zsh install-tree"  "$INS" zsh -f -c ". '$I/bbrsim_env.sh' && $Q"
# Re-sourcing another tree must switch to it (PYTHONPATH keeps both, newest first).
try "zsh re-source"     "$SRC|$I/python" zsh -f -c ". '$I/bbrsim_env.sh' && . '$R/bbrsim_env.sh' && echo \"\$BBRSIMINSTALL|\$BBRSIMLIB|\$BBRSIMINCLUDE|\$BBRSIMDATA|\${PYTHONPATH%%:*}|\${PYTHONPATH#*:}\""
try "tcsh from its dir" "$INS" tcsh -f -c "cd '$I' && source bbrsim_env.csh && $Q"
try "tcsh dir argument" "$INS" tcsh -f -c "source '$I/bbrsim_env.csh' '$I' && $Q"

T=$(mktemp -d "${TMPDIR:-/tmp}/bbrsim-test-env.XXXXXX") && T=$(cd "$T" && pwd -P) || exit 1
trap 'rm -rf "$T"' EXIT
B="$T/bogus"; W="$T/with  space/pfx"; WI="$W/share/BBRsim"   # two blanks: must not collapse
mkdir -p "$B" "$WI" && cp "$R/bbrsim_env.sh" "$R/bbrsim_env.csh" "$B/" && cp "$R/bbrsim_env.sh" "$R/bbrsim_env.csh" "$WI/" || exit 1
WINS="$WI|$W/lib|$W/include/BBRsim|$WI/data|$WI/python"
# csh: the directory argument wins over a copy in the current directory (here the repo root).
try_in "$R" "tcsh dir argument, other copy in CWD" "$INS" tcsh -f -c "source '$I/bbrsim_env.csh' '$I' && $Q"
try "tcsh source-tree"  "$SRC" tcsh -f -c "source '$R/bbrsim_env.csh' '$R' && $Q"
try "tcsh re-source"    "$SRC:$I/python" tcsh -f -c "source '$I/bbrsim_env.csh' '$I' && source '$R/bbrsim_env.csh' '$R' && $Q"
# Paths containing blanks.
try "bash blanks"       "$WINS" bash -c ". '$WI/bbrsim_env.sh' && $Q"
try "zsh blanks"        "$WINS" zsh -f -c ". '$WI/bbrsim_env.sh' && $Q"
try "tcsh blanks, from its dir" "$WINS" tcsh -f -c "cd '$WI' && source bbrsim_env.csh && $Q"
try "tcsh blanks, dir argument" "$WINS" tcsh -f -c "source '$WI/bbrsim_env.csh' '$WI' && $Q"
# Set-but-empty path variables gain no trailing ':' (an empty entry means the CWD).
E='echo "$PYTHONPATH|$LD_LIBRARY_PATH|$DYLD_LIBRARY_PATH"'
try "bash empty paths"  "$I/python|$P/lib|$P/lib" env PYTHONPATH= LD_LIBRARY_PATH= DYLD_LIBRARY_PATH= bash -c ". '$I/bbrsim_env.sh' && $E"
try "tcsh empty paths"  "$I/python|$P/lib|$P/lib" env PYTHONPATH= LD_LIBRARY_PATH= DYLD_LIBRARY_PATH= tcsh -f -c "source '$I/bbrsim_env.csh' '$I' && $E"
# A failed source (a copy outside any tree) returns 1 and changes nothing.
try "bash failed source returns 1" "rc=1" bash -c ". '$B/bbrsim_env.sh'; echo rc=\$?"
try "bash failed re-source" "$INS" bash -c ". '$I/bbrsim_env.sh'; . '$B/bbrsim_env.sh'; $Q"
try "zsh failed re-source"  "$INS" zsh -f -c ". '$I/bbrsim_env.sh'; . '$B/bbrsim_env.sh'; $Q"
try "tcsh failed re-source" "$INS" tcsh -f -c "source '$I/bbrsim_env.csh' '$I'; source '$B/bbrsim_env.csh' '$B'; $Q"
# dash cannot locate a sourced script: it uses ./bbrsim_env.sh and says so on stderr.
try_in "$I" "dash from its dir" "$INS" dash -c ". ./bbrsim_env.sh && $Q"
try_in "$R" "dash CWD fallback announced" "bbrsim_env.sh: this shell cannot locate the sourced script; using $R/bbrsim_env.sh" \
  dash -c ". '$I/bbrsim_env.sh' 2>&1 >/dev/null"
[ $bad -eq 0 ] && echo "ENV OK" || echo "ENV BAD"
exit $bad
