"""Locate the BBRsim runtime data tree (Python twin of the C++ data default).

Order: $BBRSIMDATA; else the nearest ancestor of this file holding
data/waveguides (a source checkout, an editable install, or the CMake install
<prefix>/share/BBRsim/python/bbrsim, whose ancestor share/BBRsim holds data/);
else <sys.prefix>/share/BBRsim/data. The C++ side resolves /bbr/dataDir, then
$BBRSIMDATA, then the compiled-in <install prefix>/share/BBRsim/data.

One difference: an empty $BBRSIMDATA counts as unset here, while the C++ side
takes it as given and stops with BBR011 (the env scripts never set it empty).
"""
import os
import sys


def data_dir():
    env = os.environ.get("BBRSIMDATA")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    while True:
        cand = os.path.join(here, "data")
        if os.path.isdir(os.path.join(cand, "waveguides")):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    cand = os.path.join(sys.prefix, "share", "BBRsim", "data")
    if os.path.isdir(cand):
        return cand
    raise FileNotFoundError(
        "BBRsim data directory not found: set BBRSIMDATA (source bbrsim_env.sh) "
        "to the directory that contains waveguides/ and materials/")
