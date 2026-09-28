# bbrsim_env.sh — configure a Bourne shell (sh), bash or zsh environment for BBRsim.
#
# Usage: . /path/to/bbrsim_env.sh      (repository root, or <prefix>/share/BBRsim)
#
# The script finds itself through BASH_SOURCE (bash) or ${(%):-%x} (zsh, the
# macOS default shell). BBRSIMINSTALL is recomputed on every source, so
# sourcing a second tree switches to it instead of keeping the first.
# Other sh-family shells (dash, ksh) cannot tell a sourced script its own path.
# There, cd to the script's directory first and source ./bbrsim_env.sh; the script
# prints the directory it used. Nothing is exported until the tree is recognised,
# so a failed source leaves the previous environment as it was.
#
# DYLD_LIBRARY_PATH (macOS) is searched before a binary's RPATH, and LD_LIBRARY_PATH
# (Linux) before its RUNPATH. After sourcing, every binary therefore loads libBBRsim
# from $BBRSIMLIB; in a source tree that is build/library, even for a binary built
# in build-debug/ (clang-debug preset). To use the binary's own RPATH, run it as
#   env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH <binary> <macro>

bbrsim_src=""
if [ -n "${BASH_VERSION:-}" ]; then
  bbrsim_src="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
  eval 'bbrsim_src="${(%):-%x}"'
elif [ -f ./bbrsim_env.sh ]; then
  bbrsim_src="./bbrsim_env.sh"
  echo "bbrsim_env.sh: this shell cannot locate the sourced script; using $(pwd -P)/bbrsim_env.sh" >&2
fi
if [ -z "$bbrsim_src" ]; then
  echo "ERROR: bbrsim_env.sh could not locate itself."
  echo "Please cd to the directory containing it and source it again."
  unset bbrsim_src
  return 1
fi
bbrsim_dir="$(cd "$(dirname "$bbrsim_src")" >/dev/null && pwd -P)" || bbrsim_dir=""
unset bbrsim_src

if [ -r "$bbrsim_dir/README.md" ]; then                               # source tree
  bbrsim_lib="$bbrsim_dir/build/library"
  bbrsim_inc="$bbrsim_dir/library/include"
  bbrsim_py="$bbrsim_dir/tools/python"
elif [ "$(basename "$(dirname "$bbrsim_dir")")" = "share" ]; then     # install tree
  bbrsim_top="$(dirname "$(dirname "$bbrsim_dir")")"
  bbrsim_lib="$bbrsim_top/lib"
  bbrsim_inc="$bbrsim_top/include/BBRsim"
  bbrsim_py="$bbrsim_dir/python"
  unset bbrsim_top
else
  echo "ERROR: $bbrsim_dir is neither a BBRsim source tree nor <prefix>/share/BBRsim."
  unset bbrsim_dir
  return 1
fi

export BBRSIMINSTALL="$bbrsim_dir"
export BBRSIMLIB="$bbrsim_lib"
export BBRSIMINCLUDE="$bbrsim_inc"
export LD_LIBRARY_PATH="$BBRSIMLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$BBRSIMLIB${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
export PYTHONPATH="$bbrsim_py${PYTHONPATH:+:$PYTHONPATH}"
unset bbrsim_dir bbrsim_lib bbrsim_inc bbrsim_py

# Runtime data root, read by BBRConfigManager (C++) and bbrsim.paths (Python);
# /bbr/dataDir overrides it per session.
export BBRSIMDATA="$BBRSIMINSTALL/data"
