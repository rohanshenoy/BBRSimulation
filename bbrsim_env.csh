# bbrsim_env.csh — configure a C shell (csh/tcsh) environment for BBRsim.
#
# Usage: source /path/to/bbrsim_env.csh             (interactive shell)
#        cd /path/to; source bbrsim_env.csh          (inside a script)
#        source /path/to/bbrsim_env.csh /path/to     (tcsh, inside a script)
#
# Pattern: G4CMP g4cmp_env.csh. BBRSIMINSTALL is recomputed on every source.
# The script looks for itself in the interactive command line ($_), then in the
# directory argument, then in the current directory. $_ is split at blanks, so for
# a path that contains blanks, pass the directory argument (quoted) or cd there.
# Nothing is exported until the tree is recognised, so a failed source leaves the
# previous environment as it was. csh has no "return": errors jump to the end label.
#
# DYLD_LIBRARY_PATH (macOS) is searched before a binary's RPATH, and LD_LIBRARY_PATH
# (Linux) before its RUNPATH. After sourcing, every binary therefore loads libBBRsim
# from $BBRSIMLIB; in a source tree that is build/library, even for a binary built
# in build-debug/ (clang-debug preset). To use the binary's own RPATH, run it as
#   env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH <binary> <macro>

set bbrsim_dir = ""
# Not the first command: in a sourced file tcsh's $_ still holds the previous
# command line until one command has run.
set bbrsim_args = ($_)
if ($#bbrsim_args >= 2) then
  if ("$bbrsim_args[2]" =~ */bbrsim_env.csh) then
    set bbrsim_dir = "`cd $bbrsim_args[2]:h >/dev/null && /bin/pwd -P`"
  endif
endif
unset bbrsim_args
if ("$bbrsim_dir" == "" && $#argv >= 1) then
  if (-e "$1/bbrsim_env.csh") then
    set bbrsim_dir = "`cd '$1' >/dev/null && /bin/pwd -P`"
  endif
endif
if ("$bbrsim_dir" == "") then
  if (-e bbrsim_env.csh) then
    set bbrsim_dir = "`/bin/pwd -P`"
  endif
endif
if ("$bbrsim_dir" == "") then
  echo "ERROR: bbrsim_env.csh could not locate itself."
  echo "Please cd to the directory containing it and source it again,"
  echo "or pass that directory: source <dir>/bbrsim_env.csh <dir>"
  goto bbrsim_env_end
endif

set bbrsim_up = "$bbrsim_dir:h"
if (-r "$bbrsim_dir/README.md") then                     # source tree
  set bbrsim_lib = "$bbrsim_dir/build/library"
  set bbrsim_inc = "$bbrsim_dir/library/include"
  set bbrsim_py  = "$bbrsim_dir/tools/python"
else if ("$bbrsim_up:t" == "share") then                 # install tree
  set bbrsim_top = "$bbrsim_up:h"
  set bbrsim_lib = "$bbrsim_top/lib"
  set bbrsim_inc = "$bbrsim_top/include/BBRsim"
  set bbrsim_py  = "$bbrsim_dir/python"
else
  echo "ERROR: $bbrsim_dir is neither a BBRsim source tree nor <prefix>/share/BBRsim."
  goto bbrsim_env_end
endif

setenv BBRSIMINSTALL "$bbrsim_dir"
setenv BBRSIMLIB     "$bbrsim_lib"
setenv BBRSIMINCLUDE "$bbrsim_inc"

# Prepend; an unset or empty variable gets no trailing ":" (an empty entry means the
# current directory to ld.so and to Python).
if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
if ("$LD_LIBRARY_PATH" == "") then
  setenv LD_LIBRARY_PATH "$BBRSIMLIB"
else
  setenv LD_LIBRARY_PATH "${BBRSIMLIB}:$LD_LIBRARY_PATH"
endif
if (! $?DYLD_LIBRARY_PATH) setenv DYLD_LIBRARY_PATH ""
if ("$DYLD_LIBRARY_PATH" == "") then
  setenv DYLD_LIBRARY_PATH "$BBRSIMLIB"
else
  setenv DYLD_LIBRARY_PATH "${BBRSIMLIB}:$DYLD_LIBRARY_PATH"
endif
if (! $?PYTHONPATH) setenv PYTHONPATH ""
if ("$PYTHONPATH" == "") then
  setenv PYTHONPATH "$bbrsim_py"
else
  setenv PYTHONPATH "${bbrsim_py}:$PYTHONPATH"
endif

setenv BBRSIMDATA "$BBRSIMINSTALL/data"

bbrsim_env_end:
unset bbrsim_dir bbrsim_up bbrsim_top bbrsim_lib bbrsim_inc bbrsim_py
