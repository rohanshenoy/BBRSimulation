# BBRsim

Geant4-based free-space blackbody radiation (BBR) simulation for cryogenic
experiments. Continues the work of Yen-Yung Chang (Caltech PhD, 2023, Ch. 5)
and implements the technical scope of the Golwala / Mirabolfathi NSF QIS
proposal *BBRsim: Free-Space Blackbody Radiation Simulation for
Superconducting Circuits and Cryogenic Detectors* (2025).

Superconducting-circuit qubits and sub-Kelvin phonon-mediated particle
detectors are limited by non-thermal quasiparticle populations. One cause is
free-space BBR leaking through machining gaps, cable slots, and flange-lid
joints in otherwise "sealed" cryostat chambers. BBR from a 4 K surface peaks
near λ = 1.3 mm, comparable to or larger than typical mating tolerances, so
every nominally closed chamber behaves as a leaky waveguide for long-wavelength
photons. Existing tools handle either open-volume ray-tracing or wave
propagation in gaps, not both consistently on CAD-accurate geometry.

BBRsim does both: ray-tracing in open volumes with Geant4 optical photons, and
wave propagation through gaps with pre-computed HFSS S-parameters, on the same
event-by-event footing.

## Status

**Working and validated**

- HFSS crack diffraction at normal incidence and at 45° off normal. Observed
  transmittance at 500 GHz is 50.0 % / 50.3 % (52 µm / 102 µm gaps), against
  the 50 % ideal for unpolarized light on a sub-cutoff gap.
- Frequency-keyed HFSS lookup: each photon gets the dataset nearest its
  frequency (validated on a mock five-frequency tree).
- Copper reflectance from the full complex Drude model, set by RRR and
  temperature.
- Planck thermal emitter.
- ROOT output and a Python analysis package.

**Current version (2026-09-30):** each ROOT result carries a metadata file with
its run configuration and provenance; true per-track reflection counts; batch
runs fail on a bad macro or a rejected `/bbr/` value; and Geant4 11.1 support
(see [Requirements](#requirements)).

**Implemented, unit-tested, not yet placed in a geometry:** loss-tangent dielectrics (Cirlex, Si, Ge).

**Not yet implemented:** PCB material; leakage-current analysis; calibration
geometry (BB source, mesh-TES detector); anomalous-skin-effect correction; the
upstream `G4OpBoundaryProcess` change (BBRsim wraps the stock process instead).

**Known limits:** the HFSS data cover 500 GHz only, so broadband crack results
are indicative, not quantitative. Oblique incidence is validated only on the
45° row of the HFSS grid, the only oblique row in the data.

## Repository layout

| Directory | Contents |
|---|---|
| `library/` | `libBBRsim` (CMake target `BBRsim::BBRsim`): the boundary-process wrapper, HFSS crack lookup, materials, run configuration and thermal emitter |
| `examples/testworld/` | `bbrsimTestWorld`: Cu slab with two HFSS crack volumes, Planck emitter or fixed gun |
| `examples/lightpipe/` | `bbrsimLightPipe`: 4 K → mixing-chamber light pipe, parametric or imported from `.STL` |
| `validation/` | regression fixtures, PASS/FAIL validators and the regression runner |
| `tools/` | the `bbrsim` Python package and plot scripts |
| `tests/` | C++ unit, component and regression tests (CTest) |
| `notebooks/` | physics notebooks, stored without outputs |
| `data/` | runtime data: HFSS tables, copper reference data, a sample STL |

## Quick start

```bash
git clone https://github.com/rohanshenoy/BBRSimulation.git
cd BBRSimulation
cmake --preset clang-release && cmake --build --preset clang-release && cmake --install build
. install/share/BBRsim/bbrsim_env.sh
cd examples/testworld && cmake --preset clang-release && cmake --build --preset clang-release
cd build && ./bbrsimTestWorld planck.mac
```

This writes `output/bbr.root` in `examples/testworld/build/`. To check that
everything works, set up the Python environment ([Analysis](#analysis-python))
and run `validation/Scripts/run_regression.sh` from the repository root.

The [user guide](user_guide.md) covers the test cases, every `/bbr/` command,
the output format and the analysis scripts.

## Requirements

- Geant4 11.1 or later, built with optical physics. UI and visualization
  drivers are optional. Developed and fully validated on 11.4.0 (macOS, Apple
  Clang). On 11.1.2 (Linux, GCC 13.2, the Caltech HPC) the library, the
  test-world example, all 87 C++ tests and pytest passed, as did the Planck and
  crack-transmission fixtures with their validators, all on the code before the
  2026-09-30 changes. The current sources compile against the 11.1.2 headers
  but have not yet been run there.
- The regression runner (`validation/Scripts/run_regression.sh`) assumes macOS
  and Apple Clang; on Linux, build and run the tests by hand for now.
- CMake ≥ 3.21 for the presets (3.16 for a manual configure).
- CMake ≥ 3.22 to build the tests (`BUILD_BBRSIM_TESTS`; the regression runner turns it on).
- A C++17 compiler: the same one Geant4 was built with. On macOS that is Apple
  Clang (see [Building](#building)).
- For the Python tools, a conda environment named `bbrsim`
  ([Analysis](#analysis-python)).

## Building

```bash
cmake --preset clang-release          # Apple Clang, Release, build/ -> install/
cmake --build --preset clang-release
cmake --install build
```

The top-level project builds and installs the library:

```
<prefix>/lib/libBBRsim.dylib
<prefix>/lib/cmake/BBRsim/                package config for find_package(BBRsim)
<prefix>/include/BBRsim/*.hh              library headers
<prefix>/share/BBRsim/data/               runtime data
<prefix>/share/BBRsim/bbrsim_env.{sh,csh} environment scripts, .bbrsim-version
<prefix>/share/BBRsim/python/bbrsim/      the Python package
<prefix>/bin/plot_*.py                    plot scripts
<prefix>/validation/                      copy of validation/
```

The install step is required: the executables read the HFSS data from the
installed copy by default. The prefix is `install/` in the repository. The
data path and the RPATH are fixed when the library is configured, so a
different prefix needs a reconfigure
(`cmake --preset clang-release -DCMAKE_INSTALL_PREFIX=<prefix>`), not only
`cmake --install --prefix`.

For source-tree development, build both examples with the library in one CMake
tree. This does not require an installed BBRsim library:

```bash
cmake --preset clang-release -DBUILD_BBRSIM_EXAMPLES=ON
cmake --build --preset clang-release
BBRSIMDATA="$PWD/data" build/examples/testworld/bbrsimTestWorld \
  examples/testworld/G4Macros/planck.mac
```

The source-tree examples link the current library target and use its headers.
Set `BBRSIMDATA` to the source data directory for runs that need HFSS tables;
the compiled default remains the configured install prefix. The top-level
install still installs only the library; configure an example directory
separately to install its executable. Separate example builds still use
`find_package(BBRsim)` and an installed library. In batch mode, both
executables exit nonzero if their macro fails, including a failure in a
nested macro and a `/bbr/` value the configuration rejects (for example
`/bbr/det/setCuRRR 0`): the macro stops there instead of running on with the
previous value.

`CMakePresets.json` pins `/usr/bin/clang` / `clang++` and `Release`; IDE CMake
integrations pick the preset up automatically. The manual equivalent is
`cmake -S . -B build -DCMAKE_C_COMPILER=/usr/bin/clang
-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release
-DCMAKE_INSTALL_PREFIX=install`. A `clang-debug` preset builds into
`build-debug/` and installs into `install-debug/`.

On macOS, do not let CMake pick a compiler by itself. It can select Homebrew
GCC, which compiles but fails at link against an Apple-Clang Geant4 with
undefined symbols on every API that takes `std::` types. If you see that,
check `CMAKE_CXX_COMPILER` in `build/CMakeCache.txt`.

| Option | Default | Effect |
|---|---|---|
| `WITH_GEANT4_UIVIS` | `ON` | UI and visualization; `OFF` gives a batch-only build |
| `BUILD_BBRSIM_TOOLS` | `ON` | Install the `bbrsim` Python package and the plot scripts |
| `BUILD_BBRSIM_TESTS` | `OFF` | Build the C++ tests in `tests/` and register them with CTest (`ctest --test-dir build`); the regression runner turns it on |
| `BUILD_BBRSIM_EXAMPLES` | `OFF` | Build both examples against the source-tree library target |
| `INSTALL_VALIDATION` | `ON` | Copy `validation/` into the prefix |
| `INSTALL_EXAMPLES` | `OFF` | Copy `examples/` into the prefix |

## Environment

Source the environment script before building or running the examples:

```bash
. <prefix>/share/BBRsim/bbrsim_env.sh          # bash, zsh
source <prefix>/share/BBRsim/bbrsim_env.csh    # csh, tcsh
```

Under bash, zsh and an interactive csh/tcsh the scripts find themselves, so
they work from any directory. Under dash, ksh or in a csh script, cd to the
script's directory first (a tcsh script can instead pass the directory:
`source <dir>/bbrsim_env.csh <dir>`). The repository root holds the same two
scripts, which point the variables at the source tree instead.

| Variable | Install tree | Source tree |
|---|---|---|
| `BBRSIMINSTALL` | `<prefix>/share/BBRsim` | repository root |
| `BBRSIMLIB` | `<prefix>/lib` | `build/library` |
| `BBRSIMINCLUDE` | `<prefix>/include/BBRsim` | `library/include` |
| `BBRSIMDATA` | `<prefix>/share/BBRsim/data` | `data` |

The scripts also prepend `$BBRSIMLIB` to `DYLD_LIBRARY_PATH` /
`LD_LIBRARY_PATH` and the `bbrsim` package to `PYTHONPATH`.

The data root is resolved in this order: the `/bbr/dataDir` macro command, then
`$BBRSIMDATA`, then the compiled-in `<prefix>/share/BBRsim/data`. Installed
binaries therefore find their data without any environment script.

`DYLD_LIBRARY_PATH` is searched before a binary's RPATH, so after sourcing, a
binary loads `libBBRsim` from `$BBRSIMLIB` even if it was built against another
tree. To use the binary's own RPATH, run it as
`env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH <binary> <macro>`.

## Linking your own application

```cmake
find_package(BBRsim REQUIRED)
target_link_libraries(app PRIVATE BBRsim::BBRsim)
```

Configure with `-DCMAKE_PREFIX_PATH=<prefix>` and the same compiler as the
library. Geant4 comes with it: the package config finds the Geant4 components
the library was built with, and sets `BBRsim_DATA_DIR`.
[validation/Scripts/consumer_smoke](validation/Scripts/consumer_smoke) is a
minimal working example.

## Examples

- **[testworld](examples/testworld/README.md)**, `bbrsimTestWorld`: a 4 mm Cu
  slab with two HFSS crack volumes in a 50 cm vacuum world, lit by the Planck
  emitter or a fixed photon gun.
- **[lightpipe](examples/lightpipe/README.md)**, `bbrsimLightPipe`: a 4 K →
  mixing-chamber light pipe, generated from parameters or imported from an
  ASCII `.STL`.

Each example is a standalone CMake project built against the installed
library, with its own `clang-release` preset. Its macros are copied beside the
binary, so it runs from its build directory. With no argument the binary opens
an interactive session and runs `vis.mac`. Output goes to `output/bbr.root` and
`output/bbr.metadata.json` under the directory it runs in. To start a new
simulation, copy an example directory and adapt it.

Each ROOT result has a matching `<stem>.metadata.json` containing its own category
legend, schema version, run configuration, geometry description, initial master
random-engine state, worker count, separate library and application source
fingerprints, Geant4 version, and HFSS data fingerprint. Application fingerprints
refresh when a standalone example is rebuilt, independently of the installed library.
Keep the pair together when copying results. `/analysis/setFileName results/run.root`
creates its parent directory and writes `results/run.metadata.json`; output failures
stop the run with `BBR022`. The loader still reads legacy `bbr_legend.json` files,
but new output requires its paired metadata. Master RNG state is provenance, not a
promise of identical event ordering across different worker counts.

Schema 2 adds `track_id`, `n_boundary` and `n_reflections` to both trees. The latter
two include the current boundary on crossings and the terminating boundary on
abspoints. `n_reflections` counts only reflection outcomes; `n_boundary` counts all
logged contacts. Historical `n_reflect` remains unchanged for older analysis code.
Use `(run_id, event_id, track_id)` to identify tracks within a result.
Every stock boundary status now has its own status code and event type, which
changes two things for older output: `NoRINDEX` (the photon is killed at the
boundary) is now an `absorption` event, not `other`, and the statuses older
output wrote as `Other` (`Undefined`, `Transmission`, the LUT, `Dichroic` and
coated-surface statuses) now have their own names. Compare event types across
the change with care.

The examples keep their own Geant4 user actions and delegate to the ordinary
`BBRPrimarySource`, `BBRAnalysis` and `BBRPhotonRecorder` helpers. The thermal API
now uses `ThermalSurface::InitializeSpectrum(T_K, emin_eV, emax_eV)` with read-only
`GetSpectrum()`, `GetTemperature_K()`, `GetArea()` and `GetEffArea()` accessors.
`GetBBSpecCDF` exposes `EnergyAxis()`, `PDF()` and `CDF()` as const references.
Bands that cannot be normalized (for example 10 GHz–20 THz at 0.5 mK) fail with
`BBR017` instead of emitting NaNs; the usual 4 K sampling sequence is unchanged.
Invalid box geometry or emissivity raises `BBR023`; attempting emission from
surfaces with zero total effective area raises `BBR019`.

## Validation

```bash
validation/Scripts/run_regression.sh [BUILD_DIR]     # default: build
```

The runner builds and installs the library with its C++ tests
(`BUILD_BBRSIM_TESTS=ON`), runs them (CTest), the Python tests (pytest) and
the env-script test, builds both examples against the install, runs thirteen
fixed-seed cases in parallel (seven validation fixtures, two real-tree checks
generated at run time, and four example macros)
with the Geant4 thread count pinned to 8 (`BBR_THREADS` overrides), scans their
logs for warnings, and runs the PASS/FAIL validators on their output. Last, it
runs the examples installed into a separate prefix without the env script, and
the plot scripts and the notebook's code on the fixture output. Compiler
warnings count as failures. The exit code is the number of unexpected
failures; a green run has `fail=0  xfail=1  xpass=0`. The runner also tests
batch error exits, redirected results, and the optional in-tree example build.

The one expected failure (`check_cu_serov.py`) is an open decision: the `HP_Cu`
alias (RRR 6) gives a loss 13 % below Serov's measurement. The fixtures, what
each validator checks and how to run one by hand are in
[validation/README.md](validation/README.md).

The wrapper's pass-through path is pinned by the CTest program
`testPassthrough`: stock boundary optics with and without the wrapper, fixed
seed, byte for byte, also with WLS active. Run the runner before merging any
change to the wrapper.

## Analysis (Python)

One-time setup, from the repository root:

```bash
conda create -n bbrsim python=3.11 numpy scipy matplotlib pandas
conda run -n bbrsim pip install -e "tools/python[test]"     # the bbrsim package, uproot and pytest
```

Run every script as `conda run -n bbrsim python <script>`, not `python3`. The
editable install ties the `bbrsim` environment to this checkout; for another
checkout, source its env script or rerun `pip install -e` there. The package's
tests are `conda run -n bbrsim python -m pytest -q -p no:cacheprovider tools/python/tests`,
run from the repository root (not from `tools/python`, where the source package
would shadow the one on `PYTHONPATH`).

```python
from bbrsim.io import load
crossings, abspoints = load("output/bbr.root")   # decoded DataFrames
```

The `bbrsim` package (`io`, `physics`, `select`, `hfss`, `paths`) holds the
ROOT loader and the only Python copy of the physics formulas; every validator
and plot script uses it. Plot scripts, installed to `<prefix>/bin`:

- `plot_cu_reflectance.py`: Cu reflectance vs frequency and temperature
- `plot_crack_angular.py`: crack exit angles vs the HFSS far field
- `plot_crack_frequency.py`: the frequency-keyed lookup
- `plot_test_output.py`: overview of one output file

[notebooks/copper_reflectance.ipynb](notebooks/copper_reflectance.ipynb) walks
the copper Drude model end to end. A Jupyter kernel needs
`conda run -n bbrsim pip install ipykernel` once.

## Physics

**Open volumes: ray tracing.** Photons are Geant4 optical photons. At metal
surfaces, BBRsim replaces Fresnel optics, which is unreliable for cryogenic
metals in this band, with a tabulated `REFLECTIVITY` stored on the material.
For copper the table comes from the full complex Drude model, parameterized by
RRR and temperature. Photons reflect specularly or are absorbed.

**Narrow gaps: wave propagation.** A gap is a volume made of the flag material
`vacuum_wg`, named after its HFSS dataset. When a photon enters one, BBRsim
samples transmission or reflection, the exit direction, polarization and exit
position from ANSYS HFSS simulations of that gap. It picks the dataset nearest
the photon's frequency; there is no interpolation between frequencies.

**Both are handled by one class.** `BBSimOpBoundaryProcess` wraps the stock
Geant4 `G4OpBoundaryProcess`. It handles photons entering `vacuum_wg` volumes
or `REFLECTIVITY` materials and passes every other boundary to the unmodified
stock process.

**Emission.** A box surface emits photons with energies drawn from the Planck
photon-number spectrum (10 GHz–20 THz), uniformly in θ over the outward
hemisphere, weighted by area and emissivity.

**Dielectrics.** Cirlex, Si and Ge have a constant refractive index and an
absorption length set by their loss tangent; the stock Geant4 processes handle
them.

The copper model and its validation against Serov (2016) are derived in the
notebook. How each physics claim is tested is in
[validation/README.md](validation/README.md).

## Versioning

Every build writes `git describe --always --dirty` to `.bbrsim-version`,
installed to `<prefix>/share/BBRsim/`. Outside a git checkout the project
version is used instead. Releases are annotated tags of the form
`bbrsim-VXX-YY-ZZ`, each with a line in [ChangeHistory](ChangeHistory).

## Upgrading an older checkout

Checkouts from before `bbrsim-V00-01-00` used a different layout: delete the old `build/`, rebuild and install as above, and run `./bbrsimTestWorld` from `examples/testworld/build/`; it replaces the executable the old layout built in `build/`.

## Software license

To be determined (PI decision; required before public release). Portions
derived from Geant4 example code remain under the Geant4 Software License.

## History

BBRsim began as a modification of the Geant4 OpNovice2 example.

## Project plan and team

From the NSF proposal, period of performance 2026Q4–2029Q3:

- **O1**: complete BBRsim development (this repository).
- **O2**: validate against a temperature-controlled OFHC-Cu BB source (4–20 K)
  and purpose-built mesh-TES photon detectors, at TAMU's DR and on the SuperCDMS
  SNOLAB Pathfinder tower at SLAC.
- **O3**: release BBRsim as a Geant4 module; upstream the
  `G4OpBoundaryProcess` `REFLECTIVITY` change through the Geant4 EM Physics
  Working Group.

Team: Golwala (PI, Caltech), Mirabolfathi (co-PI, TAMU), Shenoy (student lead,
Caltech), Xiong (Brinson postdoc, Caltech), Kurinsky / Partridge (SLAC
SuperCDMS, unfunded collaborators), Chang (BBRsim originator, unfunded advisor).

## References

- Chang, Y.-Y. (2023). *SuperCDMS HVeV Run 2 Low-Mass Dark Matter Search,
  Highly Multiplexed Phonon-mediated Particle Detector with Kinetic Inductance
  Detector, and the Blackbody Radiation in Cryogenic Experiments*. PhD thesis,
  Caltech. Chapter 5 is the primary physics reference.
- Golwala & Mirabolfathi, *BBRsim* NSF QIS proposal (2025).
- Agostinelli et al., "GEANT4 — a simulation toolkit," *NIM A* **506**, 250 (2003).
- Griffiths, D.J. (2017). *Introduction to Electrodynamics*, 4th ed. §9.4.
- Serov, Y.L. et al. (2016). *IEEE Trans. Microwave Theory Tech.* **64**(11), 3828.
- Poole et al., CADMesh (2nd ver.): CAD → Geant4 geometry import.
