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
propagation in gaps — not both consistently on CAD-accurate geometry.

BBRsim does both: ray-trace in open volumes via Geant4 optical photons, and
wave propagation through gaps via pre-computed HFSS S-parameters — on the same
event-by-event footing, so BBR backgrounds can be simulated rather than
debugged after the fact.

The repository follows the G4CMP layout:

| Directory | Contents |
|---|---|
| `library/` | `libBBRsim` (CMake target `BBRsim::BBRsim`): the boundary-process wrapper, HFSS crack lookup, materials, config manager and thermal emitter |
| `examples/testworld/` | `bbrsimTestWorld`: Cu slab with two HFSS crack volumes, Planck emitter or fixed gun |
| `examples/lightpipe/` | `bbrsimLightPipe`: 4 K → mixing-chamber light pipe, parametric or CAD (`.STL`) |
| `validation/` | `Validation_*.mac` fixtures, `check_*.py` validators, the regression runner |
| `tools/` | the `bbrsim` Python package (`tools/python/bbrsim`) and the `plot_*.py` scripts |
| `notebooks/` | physics notebooks, stored without outputs |
| `data/` | runtime data (HFSS tables, material reference data, a sample STL), installed to `share/BBRsim/data` |

## Software License

To be determined (PI decision; required before public release). Portions
derived from Geant4 example code remain under the Geant4 Software License.

## Downloading

```bash
git clone https://github.com/rohanshenoy/BBRSimulation.git
cd BBRSimulation
```

Release tags have the form `bbrsim-VXX-YY-ZZ` and are listed in
[ChangeHistory](ChangeHistory); check one out with `git checkout
bbrsim-VXX-YY-ZZ`. Then build as described under
[Building the Package](#building-the-package).

## User Environment

Requirements: Geant4 11.x built with optical physics (developed against
11.4.0), CMake ≥ 3.21 for the presets used below (3.16 for a manual
configure), a C++17 compiler (Apple Clang on macOS, see below). The Python
tools need the `bbrsim` conda environment (see
[Analysis (Python)](#analysis-python)).

After installing, source the BBRsim environment script before building or
running the examples:

```bash
. <prefix>/share/BBRsim/bbrsim_env.sh          # bash, zsh (dash, ksh: see below)
source <prefix>/share/BBRsim/bbrsim_env.csh    # csh, tcsh
```

`<prefix>` is the install prefix, `install/` in the repository by default.
Under bash and zsh, and in an interactive csh/tcsh, the scripts locate
themselves, so they can be sourced from any directory. Other sh-family shells
(dash, ksh) and csh scripts cannot: cd to the script's directory and source it
from there (a tcsh script can instead pass the directory,
`source <dir>/bbrsim_env.csh <dir>`). Sourcing a second tree switches every
variable to it. The repository root holds the same two scripts; sourcing those
points the variables at the source tree instead.

| Variable | Install tree | Source tree |
|---|---|---|
| `BBRSIMINSTALL` | `<prefix>/share/BBRsim` (the script's directory) | repository root |
| `BBRSIMLIB` | `<prefix>/lib` | `build/library` |
| `BBRSIMINCLUDE` | `<prefix>/include/BBRsim` | `library/include` |
| `BBRSIMDATA` | `<prefix>/share/BBRsim/data` | `data` |

The scripts also prepend `$BBRSIMLIB` to `DYLD_LIBRARY_PATH` and
`LD_LIBRARY_PATH`, and the `bbrsim` package directory
(`<prefix>/share/BBRsim/python`, or `tools/python` in the source tree) to
`PYTHONPATH`.

`BBRSIMDATA` is the runtime data root (it must contain `waveguides/`) and is the
environment counterpart of the `/bbr/dataDir` macro command. Precedence:
`/bbr/dataDir` (in a macro, before `/run/initialize`) > `$BBRSIMDATA` > the
compiled-in default `<prefix>/share/BBRsim/data`. Installed binaries therefore
find their data without any environment script.

`DYLD_LIBRARY_PATH` (macOS) and `LD_LIBRARY_PATH` (Linux) are searched before a
binary's RPATH. After sourcing, every binary loads `libBBRsim` from
`$BBRSIMLIB`; with the source-tree script that is `build/library`, even for a
binary built against another tree. To use a binary's own RPATH, run it as
`env -u DYLD_LIBRARY_PATH -u LD_LIBRARY_PATH <binary> <macro>`.

## Building the Package

```bash
cmake --preset clang-release          # Apple Clang, Release, build/ -> install/
cmake --build --preset clang-release
cmake --install build
```

The top-level project builds and installs the library only:

```
<prefix>/lib/libBBRsim.dylib
<prefix>/lib/cmake/BBRsim/                package config for find_package(BBRsim)
<prefix>/include/BBRsim/*.hh              the 12 library headers
<prefix>/share/BBRsim/data/               runtime data
<prefix>/share/BBRsim/bbrsim_env.{sh,csh} .bbrsim-version
<prefix>/share/BBRsim/python/bbrsim/      the Python package
<prefix>/bin/plot_*.py                    plot scripts
<prefix>/validation/                      copy of validation/
```

Install goes to `./install` by default. The data default and the RPATH are
baked in at configure time, so changing the prefix requires reconfiguring
(`cmake --preset clang-release -DCMAKE_INSTALL_PREFIX=<prefix>`); `cmake
--install build --prefix <other>` alone installs a library whose data default
and RPATH still name the configured prefix. The install step is required: with
`BBRSIMDATA` unset and no
`/bbr/dataDir`, the executables read the HFSS data from the installed copy, so
a crack macro run before installing aborts with `BBR011`.

`CMakePresets.json` pins `/usr/bin/clang`/`clang++` and `Release`; IDE CMake
integrations (VSCode CMake Tools) pick the preset up automatically instead of
offering their own compiler kits. The manual equivalent is
`cmake -S . -B build -DCMAKE_C_COMPILER=/usr/bin/clang
-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release
-DCMAKE_INSTALL_PREFIX=install`; a `clang-debug` preset builds into
`build-debug/` and installs into `install-debug/`.

On macOS, configure with Apple Clang explicitly as shown. A bare `cmake ..` can
pick up Homebrew GCC (libstdc++), which compiles but fails at link against a
Geant4 built with Apple Clang (libc++) on every API whose signature contains
`std::` types. If you see undefined Geant4 symbols at link time, check
`CMAKE_CXX_COMPILER` in `build/CMakeCache.txt`.

Options (pass as `-D<option>=ON|OFF` to the configure step):

| Option | Default | Effect |
|---|---|---|
| `WITH_GEANT4_UIVIS` | `ON` | Build with Geant4 UI and visualization drivers; `OFF` gives a batch-only build |
| `BUILD_BBRSIM_TOOLS` | `ON` | Install the `bbrsim` Python package and the plot scripts |
| `INSTALL_VALIDATION` | `ON` | Copy `validation/` into the prefix |
| `INSTALL_EXAMPLES` | `OFF` | Copy `examples/` into the prefix |

Batch-only (no UI/visualization):

```bash
cmake --preset clang-release -DWITH_GEANT4_UIVIS=OFF
cmake --build --preset clang-release
cmake --install build
```

## Linking user applications

A CMake project finds the installed library with `find_package` and links the
imported target. Geant4 comes with it: the package config replays the Geant4
components the library was built with.

```cmake
find_package(BBRsim REQUIRED)
target_link_libraries(app PRIVATE BBRsim::BBRsim)
```

Configure with `-DCMAKE_PREFIX_PATH=<prefix>` (and, on macOS, the same Apple
Clang as the library). The package config also sets `BBRsim_DATA_DIR`.
[validation/Scripts/consumer_smoke](validation/Scripts/consumer_smoke) is the
minimal example: one CMakeLists.txt and one source file that includes
`BBRConfigManager.hh` and `BBSimPhysics.hh` and prints the compiled-in data
default. The regression runner builds and runs it on every pass.

## Application Examples

- **testworld** ([examples/testworld](examples/testworld/README.md)) —
  `bbrsimTestWorld`: a 4 mm Cu slab with two HFSS crack volumes in a 50 cm
  vacuum world, lit by the Planck thermal emitter or a fixed photon gun. The
  Planck, Cu reflectance and crack runs below use it.
- **lightpipe** ([examples/lightpipe](examples/lightpipe/README.md)) —
  `bbrsimLightPipe`: a 4 K → mixing-chamber light pipe, generated from
  parameters or imported from an ASCII `.STL`.

Each example is a standalone CMake project built against the installed
library, with its own `CMakePresets.json` (`clang-release`: build dir `build/`;
library search path: the prefix of a sourced env script, else `../../install`).
The macros in `G4Macros/` are copied beside the binary, so it runs from its
build directory:

```bash
. install/share/BBRsim/bbrsim_env.sh
cd examples/testworld && cmake --preset clang-release && cmake --build --preset clang-release
cd build && ./bbrsimTestWorld planck.mac
```

With no argument the binary starts an interactive session and runs `vis.mac`
(needs a Geant4 build with UI and visualization drivers). Output goes to
`output/bbr.root` and `output/bbr_legend.json` under the directory the binary
runs in. An example directory can be copied elsewhere and adapted; source the
env script (or pass `-DCMAKE_PREFIX_PATH=<prefix>`) so it finds the library.

## Validation

```bash
validation/Scripts/run_regression.sh [BUILD_DIR]     # default: build
```

The runner configures, builds and installs the library (prefix `BBR_PREFIX`,
default `install/`), sources the installed env script, builds both examples
against it, runs the drift guards and the consumer smoke test, runs the eight
fixtures in parallel (each in `BUILD_DIR/regression/<case>/`), scans their logs
for `GeomNav`, `G4Exception`, `BBR0xx` and `LP002` messages, and then runs the
validators on their output. It prints one PASS/FAIL line per check; the exit
code is the number of unexpected failures, and a compiler warning counts as a
failure (the builds are incremental, so only files compiled in that run are
seen: use a fresh `BUILD_DIR` for a full warnings audit). Overrides:
`BBR_PREFIX`, `BBR_PYTHON`, `BBR_JOBS`, `BBR_XFAIL`. A green run ends with the
line `pass=41  fail=0  xfail=1  xpass=0` (two spaces between fields, then the
output directory).

| Fixture (`validation/G4Macros/`) | Executable | Validators | Passes when |
|---|---|---|---|
| `Validation_Reflectance.mac` | `bbrsimTestWorld` | `check_reflectance.py`, `check_term_status.py`, `check_no_photons_in_metal.py` | absorbed count within 5 σ (Poisson) of full-Drude theory |
| `Validation_Planck.mac` | `bbrsimTestWorld` | `check_planck_spectrum.py --temp 4`, `check_nreflect.py`, `check_angle_distribution.py`, `check_term_status.py`, `check_no_photons_in_metal.py` | spectral peak within [0.65, 1.35] of the photon-number peak; n_reflect = 1 modal and falling; incidence angles uniform in θ (KS) |
| `Validation_CrackWall.mac` | `bbrsimTestWorld` | `check_crack_wall_reflection.py`, `check_no_photons_in_metal.py`, `check_term_status.py` | crack-wall reflections flip p_z only; no photon inside Cu |
| `Validation_WorldExit.mac` | `bbrsimTestWorld` | `check_term_status.py` | every world exit is labelled `WorldExit` |
| `Validation_CrackTransmit.mac` | `bbrsimTestWorld` | `check_crack_transmittance.py`, `check_no_photons_in_metal.py`, `check_term_status.py` | T_obs within 3 σ of 0.50; no tangential exits |
| `Validation_CrackOblique.mac` | `bbrsimTestWorld` | `check_crack_oblique.py`, plus the two invariants on each of its 16 files | all 138 oblique-incidence checks |
| `Validation_CrackFrequency.mac` | `bbrsimTestWorld` | `check_crack_frequency.py`, plus the two invariants on each of its 16 files | all 89 frequency-lookup checks; one `BBR008` clamp warning per side |
| `Validation_LightPipe.mac` | `bbrsimLightPipe` | `check_no_photons_in_metal.py`, `check_term_status.py` | no photon inside the wall; every termination labelled |
| — (no ROOT input) | — | `check_physics.py`, `check_cu_serov.py` | Python formulas match documented reference values; Serov points within ±10 % |
| — (run by hand) | — | `check_crack_ratio.py`, `check_cu_absorptance.py` | crack2/crack1 entry ratio vs aperture ratio; Planck-weighted Cu absorptance |

That is 8 fixtures and 14 validators: the runner calls 12, and
`check_crack_ratio.py` and `check_cu_absorptance.py` are run by hand on a
suitable output. `check_cu_serov.py` is an expected failure (XFAIL): the
`HP_Cu` alias (RRR 6) gives a loss 13 % below Serov's measurement under the
full Drude model, and the choice between RRR 5 and a wider tolerance is still
open. The runner also calls `validation/Scripts/drift_guards.sh`, which checks
that the CMake file lists match the files on disk, that the library names no
example class, that every fixture is run, and that tracked notebooks carry no
outputs. Details, the mock HFSS tree and the manual validators:
[validation/README.md](validation/README.md).

**Wrapper regression.** `BBSimOpBoundaryProcess` falls through to the stock
`G4OpBoundaryProcess` for any geometry without `vacuum_wg` or `REFLECTIVITY`
materials. The historical fixed-seed pass-through macro (`verify_wrapper.mac`)
was retired with the old example geometry. Re-establishing an automated
pass-through regression, on a purpose-built minimal geometry that exercises
stock boundary optics, is open work. Until then, run the regression runner
before merging any change to the wrapper.

## Analysis (Python)

The `bbrsim` package (`tools/python/bbrsim`: `io`, `physics`, `select`,
`hfss`, `paths`) is the single source of truth for reading BBRsim output and for
the Drude / Planck / Hagen-Rubens formulas; every `check_*` and `plot_*` script
reads through it. One-time setup, from the repository root:

```bash
conda create -n bbrsim python=3.11 numpy scipy matplotlib pandas
conda run -n bbrsim pip install -e tools/python     # bbrsim + uproot, the ROOT reader
```

Both steps are needed: `bbrsim.io` imports uproot, which only the editable
install brings in. Sourcing either env script also puts `bbrsim` on
`PYTHONPATH`; the editable install makes it importable without the env script
(a Jupyter kernel, say). The editable install points the shared `bbrsim` env at
this checkout: for another clone or an install, source that tree's env script
(`PYTHONPATH` takes precedence) or rerun `pip install -e` there.

Run every script as `conda run -n bbrsim python <script>`, not `python3`.
`bbrsim.paths.data_dir()` is the Python twin of the C++ data default:
`$BBRSIMDATA`, else the nearest `data/` holding `waveguides/` above the package,
else `<sys.prefix>/share/BBRsim/data`. One difference: an empty `BBRSIMDATA`
counts as unset in Python, while the C++ side takes it as given and stops with
`BBR011`.

```python
from bbrsim.io import load
crossings, abspoints = load("output/bbr.root")   # decoded DataFrames
```

Plot scripts live in `tools/` and are installed to `<prefix>/bin`:
`plot_cu_reflectance.py` (Cu reflectance vs frequency and temperature),
`plot_crack_angular.py` (crack exit angles vs the HFSS far field),
`plot_crack_frequency.py` (the frequency-keyed lookup) and
`plot_test_output.py` (overview of a crossings file). `plot_cu_reflectance.py`
writes to the current directory (override with `--out`);
`plot_test_output.py` and `plot_crack_angular.py` write next to their input
file; `plot_crack_frequency.py` writes into its output directory (or `--out`).
`notebooks/copper_reflectance.ipynb` walks the copper Drude
model end to end against the same package (a kernel needs
`conda run -n bbrsim pip install ipykernel` once).

## Versioning

Every build writes `git describe --always --dirty` to `.bbrsim-version`, which
is installed as `<prefix>/share/BBRsim/.bbrsim-version`. When the source
directory is not the top of its own git work tree (an unpacked tarball, or a
copy vendored inside another repository), the `BBRSIM_VERSION` cache variable,
default the project version, is used instead.
Release tags are annotated, because `git describe` ignores lightweight tags,
and have the form `bbrsim-VXX-YY-ZZ`; each gets a line in
[ChangeHistory](ChangeHistory).

## Physics

### Status (September 2026)

Core physics is operational and validated. The HFSS diffraction path, the Cu
reflectance model, the Planck thermal emitter, and the ROOT output/analysis
layer are implemented and tested. The loss-tangent dielectrics (Cirlex, Si, Ge)
are implemented but not yet placed in any geometry.

As of 2026-08-26 the `fix/core-hardening` and `light-pipe-example` branches are
**merged into `main`**, so everything below ships from a single branch: the core
simulation, the navigator/cache hardening, and the light-pipe example with CAD
(`.STL`) import. In September 2026 the repository was reorganized into the
library / examples / validation / tools layout above, with behaviour unchanged.

#### Recent correctness hardening

- HFSS diffraction now relocates Geant4's navigator to a point just inside the
  crack before applying a non-local exit state. This prevents stale safety and
  touchable state after the in-volume transport shortcut.
- The shared crack/HFSS dataset cache is protected during lazy initialization
  and lookup, making concurrent worker access safe in multithreaded runs. It is
  shared mutable state with synchronized lazy initialization, not an immutable
  singleton.
- CADMesh's optional reverse-coordinate flag is explicitly initialized, so CAD
  light-pipe construction does not depend on indeterminate state.
- (September 2026 audit) The tabulated-reflectance handler now takes the
  reflecting normal from the navigator's exit normal, as stock Geant4 does. It
  previously used the entered solid's `SurfaceNormal()`, which for a photon
  leaving a crack through its side wall into the Cu slab returned the slab's
  x-face and sent the photon through solid copper (`Validation_CrackWall.mac`).
- `abspoints.term_status` is now `WorldExit` for photons killed at the world
  boundary; it used to echo the boundary process's stale per-thread status,
  including `BBRAbsorb` for photons that were never absorbed
  (`Validation_WorldExit.mac`).
- The Planck emitter box is configurable (`/bbr/thermal/emitterCenter`,
  `/bbr/thermal/emitterSize`); `lightpipe.mac` sizes it to the bore, so photons
  are no longer created inside the light-pipe wall.
- Raw HFSS power ratios above 1 are capped at load time and exit-face-plane
  far-field directions get zero sampling weight; crack transmittance at normal
  incidence is 50.0% / 50.3% against the 50% ideal
  (`Validation_CrackTransmit.mac`).

These changes were smoke-tested with a 10,000-event fixed-gun run using 15
workers; no geometry-navigation warnings, boundary-process errors, or stuck
tracks were observed. There is no registered `ctest` suite: correctness is
checked by the `validation/check_*.py` PASS/FAIL validators, which the
regression runner drives (see [Validation](#validation)).

#### Working

- **HFSS diffraction** — `BBRHFSSData` loads far-field + waveguide CSVs for one
  frequency; `BBRCrackLibrary` discovers each crack's frequency grid from the
  dataset directory names and lazy-loads the grid point nearest the photon;
  `BBSimOpBoundaryProcess` intercepts photons entering `vacuum_wg` crack
  volumes and routes them through the HFSS lookup.
  Validated at 500 GHz normal incidence: observed transmittance
  50.0 ± 0.25% (52 µm gap) / 50.3 ± 0.25% (102 µm gap) over 40k photons each,
  against the exact 50% unpolarized ideal — only the TEM component transmits
  through a sub-cutoff gap. Raw HFSS power ratios above 1 (a port-normalization
  artefact, 1.0545 for the 52 µm gap) are capped at 1 when the tables load;
  before that cap the observed values were 51.9% / 50.5%.

- **Cu reflectance** — Full complex Drude model (σ(ω) = σ_DC/(1−iωτ) in
  ε̃ = 1 + iσ/(ε₀ω); Griffiths §9.4 generalized) parameterized by RRR and
  temperature. `BBRMaterials::GetCopper(RRR, T_K)` builds a 24-point log-spaced
  REFLECTIVITY table from 10 GHz to 20 THz. Three named grades available; users
  may also supply any integer RRR directly.

- **Planck thermal emitter** — `ThermalSurface` + `GetBBSpecCDF`: box-surface
  emitter sampling the Planck photon-number spectrum (10 GHz–20 THz), with
  per-surface emissivity weighting. Default mode of the examples' primary
  generator (`TestWorldPrimaryGeneratorAction` and its light-pipe copy);
  temperature set at runtime via `/bbr/thermal/setT`.

- **Test geometry** — `TestWorldDetectorConstruction` (`examples/testworld`):
  50 cm world, 4 mm Cu slab, two `vacuum_wg` crack daughters at z = 0 and
  z = 3 mm. Cu material is configurable via `/bbr/det/` before
  `/run/initialize`; gun and emitter settings via `/bbr/gun/` and
  `/bbr/thermal/` at any time.

- **Loss-tangent dielectrics** — `BBRMaterials::GetCirlex()`,
  `GetSiliconCrystal()`, `GetGermaniumCrystal()`: flat `RINDEX` plus
  `ABSLENGTH = c/(2πν·n·tanδ)` on the same 10 GHz–20 THz grid, with bulk
  absorption via stock `G4OpAbsorption` and Fresnel via the stock boundary
  process. Compile- and link-verified only — not yet placed in any geometry
  and not covered by a `check_*` validator.

- **ROOT output + analysis layer** — each example's run and stepping actions
  (`TestWorldRunAction` / `TestWorldSteppingAction`, and the light-pipe copies)
  write `output/bbr.root` (`crossings` + `abspoints` ntuples) and
  `output/bbr_legend.json` via `G4AnalysisManager` with ntuple merging under
  multithreading. `notebooks/copper_reflectance.ipynb` walks the copper Drude
  model end to end against the same package. `tools/python/bbrsim/` (`io.py`,
  `physics.py`, `select.py`, `hfss.py`, `paths.py`) is the single source of
  truth for loading and for the Drude / Planck / Hagen-Rubens formulas; every
  `check_*` and `plot_*` script reads through it.

#### Not yet implemented

- Patched `G4OpBoundaryProcess` with `REFLECTIVITY` on `G4MaterialPropertiesTable`
  (upstream Geant4 PR target — BBRsim currently *wraps* the stock process rather
  than patching it)
- Leakage-current post-processing (the `abspoints` ntuple that feeds it exists;
  the loss-tangent-weighting analysis is Phase B2)
- PCB material; geometry placement and a transmission validator for the
  existing Si/Ge/Cirlex getters
- BBR calibration geometry (BB source + mesh-TES detector models)
- Anomalous-skin-effect correction to the classical Drude boundary optics

#### Known scope limits

- **HFSS data is single-frequency; the lookup is not.** The dataset is chosen
  per photon by nearest frequency in log space, and the choice is recorded in
  the `hfss_freq_GHz` output column. Only a 500 GHz dataset exists per crack,
  so in practice every photon still gets the 500 GHz tables and broadband crack
  results stay indicative rather than quantitative. Dropping in real exports
  means adding `<id>_<freq>GHz_Ephi=N` directories; there is no interpolation
  between grid points. The mechanism is validated against a mock five-frequency
  tree (`Validation_CrackFrequency.mac`).
- **Oblique incidence is validated at 45° only.** `Validation_CrackOblique.mac`
  and `check_crack_oblique.py` verify the azimuth fold/unfold on the θ = 135°
  row of the HFSS grid (the only oblique incidence in the data): momentum along
  the plate's long axis keeps its sign, the polarization filter behaves (E
  across the gap transmits, in-plane E is cut off), ±y and ±z tilts mirror.
  The z-mirror is unobservable on a parallel plate, and no HFSS run outside
  the [0°, 90°] wedge exists yet.

### Planck emitter and crack diffraction

From `examples/testworld/build/`, after sourcing the env script:

```bash
./bbrsimTestWorld planck.mac          # 10 000 thermal photons at 4 K
./bbrsimTestWorld planck_5M.mac       # 5M-event production run at 4 K
./bbrsimTestWorld planck_10K.mac      # 1M events at 10 K
```

Photons entering the `vacuum_wg` cracks are routed through the HFSS lookup;
expected transmittance at 500 GHz normal incidence is 50% per crack (the
unpolarized TEM-only ideal; observed 50.0% / 50.3%, see
`Validation_CrackTransmit.mac` and `check_crack_transmittance.py`). Output is a
ROOT file `output/bbr.root` — two ntuples, `crossings` (one row per
optical-photon boundary crossing) and `abspoints` (one row per photon
termination) — plus a `output/bbr_legend.json` sidecar that maps the integer
code columns (status / event_type / volume / material) back to names. Runs are
multithreaded; G4Analysis merges the per-thread ntuples. Read it in Python via
the shared loader `tools/python/bbrsim/io.py`, which decodes the codes to the
legacy column names:

```bash
conda run -n bbrsim python <repo>/validation/check_planck_spectrum.py output/bbr.root --temp 4
```

Every `check_*` / `plot_*` script reads BBRsim output through the `bbrsim`
package; none of them read a BBRsim-produced CSV. The only remaining
`pandas.read_csv` calls load external reference data (HFSS far-field tables and
the Palik / Serov / Geant4-IR copper comparison sets).

### Cu reflectance

```bash
./bbrsimTestWorld reflectance.mac          # OFHC_Cu (RRR=100), 500 GHz, 10 000 events
./bbrsimTestWorld reflectance_OF_Cu.mac    # OF_Cu   (RRR=3),   500 GHz, 2 000 events
./bbrsimTestWorld reflectance_HP_Cu.mac    # HP_Cu   (RRR=6),   500 GHz, 2 000 events
```

Output: `output/bbr.root`, checked by `check_reflectance.py` below. Each worker
thread also keeps a running tally, printed as `G4WTn > [BBR] reflectance
mat=... N=... A_obs=... R_theory=...` every 1000 Cu hits on that thread, so a
10 000-event run spread over many threads may print none; add
`/run/numberOfThreads 1` before `/run/initialize` for a running tally.
At 4 K OFHC Cu sits on the relaxation plateau (D ≈ 4.9×10⁻⁵ at 500 GHz).
Compare against Drude theory:

```bash
./bbrsimTestWorld reflectance.mac
conda run -n bbrsim python <repo>/validation/check_reflectance.py
```

Plot reflectance vs frequency across Cu grades and temperatures:

```bash
conda run -n bbrsim python <repo>/tools/plot_cu_reflectance.py
# output: cu_reflectance_plots.png in the current directory (override with --out)
```

### Setting Cu material at runtime

In any mac file, before `/run/initialize`:

```mac
# Named alias (sets RRR and temperature to 4 K)
/bbr/det/setCuMaterial OFHC_Cu

# Direct RRR input (any integer >= 1)
/bbr/det/setCuRRR 300

# Warm shield layer: set temperature first, then RRR
/bbr/det/setCuStageT 40 K
/bbr/det/setCuRRR 50
```

Valid aliases: `OFHC_Cu` (RRR=100), `OF_Cu` (RRR=3), `HP_Cu` (RRR=6).
σ_DC is derived automatically as RRR × 5.96×10⁷ S/m.

### Light-pipe geometry

The light-pipe example, `bbrsimLightPipe`, models a 4 K → mixing-chamber light
pipe. It links the same library (physics list, materials, HFSS lookup, emitter)
and carries its own copies of the test world's primary-generator, stepping, run
and action-initialization classes; only the detector construction and its
messenger differ.

```bash
cd examples/lightpipe && cmake --preset clang-release && cmake --build --preset clang-release
cd build && ./bbrsimLightPipe lightpipe.mac
```

`lightpipe.mac` places the Planck emitter just upstream of the warm aperture
and inside the 5 mm bore (`/bbr/thermal/emitterCenter -51 0 0 mm`,
`/bbr/thermal/emitterSize 1 7 7 mm`). With the test-world default emitter (a
1×20×20 mm patch centred at x = −50 mm) the emitting face sits inside the tube
wall and photons are created in the copper; `BuildParametric` warns (`LP002`)
if the configured emitter reaches into the wall.

Two build modes, selected by `/bbr/lightpipe/mode`:

- **`parametric`** (default) — a tube generated from bore radius, length, and
  wall thickness.
- **`cad`** — an ASCII `.STL` imported through the bundled header-only CADMesh
  (`library/include/CADMesh.hh`). The built-in reader is ASCII-only; binary STL
  needs assimp. A sample mesh ships at `data/cad/box_sample.stl`.

| Command | Argument | Description |
|---|---|---|
| `/bbr/lightpipe/mode` | `parametric \| cad` | Build mode |
| `/bbr/lightpipe/bore` | length + unit | Inner bore radius (aperture) |
| `/bbr/lightpipe/length` | length + unit | Tube length along +x |
| `/bbr/lightpipe/wallThickness` | length + unit | Wall thickness |
| `/bbr/lightpipe/wallMaterial` | `Cu \| reflector` | Wall optical material |
| `/bbr/lightpipe/stlPath` | path | ASCII `.STL` to load (`cad` mode) |

All `/bbr/lightpipe/` commands are **`PreInit` only** and are not broadcast to
worker threads — issue them before `/run/initialize`. There is no runtime
geometry reinitialization; changing a parameter after initialization requires a
new session. The differential 4 K → MXC measurement this geometry is meant to
support is not implemented yet.

## Migration from the pre-reorg layout

After pulling, delete stale `build/` and `build-audit/` trees, then rebuild and install as above; `./BBRSim planck.mac` becomes `./bbrsimTestWorld planck.mac`, run from `examples/testworld/build/`.
The light-pipe executable is now `bbrsimLightPipe`, the user macros live in
`examples/*/G4Macros/`, the regression macros are
`validation/G4Macros/Validation_*.mac`, the validators and the runner are under
`validation/`, and the plot scripts and the Python package are under `tools/`.

## History

BBRsim began as a modification of the Geant4 OpNovice2 example.

## Project plan and team

From the NSF proposal, period of performance 2026Q4–2029Q3:

- **O1** — Complete BBRsim development (this repository).
- **O2** — Validate against a temperature-controlled OFHC-Cu BB source (4–20 K)
  plus purpose-built mesh-TES photon detectors at TAMU's DR and the SuperCDMS
  SNOLAB Pathfinder tower at SLAC.
- **O3** — Release BBRsim as a Geant4 module; upstream the `G4OpBoundaryProcess`
  `REFLECTIVITY` change via the Geant4 EM Physics Working Group.

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
- Poole et al., CADMesh (2nd ver.) — CAD → Geant4 geometry import.
