# BBRsim

Geant4-based free-space blackbody radiation (BBR) simulation for cryogenic
experiments. Continues the work of Yen-Yung Chang (Caltech PhD, 2023, Ch. 5)
and implements the technical scope of the Golwala / Mirabolfathi NSF QIS
proposal *BBRsim: Free-Space Blackbody Radiation Simulation for
Superconducting Circuits and Cryogenic Detectors* (2025).

## Why this exists

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

## Status (September 2026)

Core physics is operational and validated. The HFSS diffraction path, the Cu
reflectance model, the Planck thermal emitter, and the ROOT output/analysis
layer are implemented and tested. The loss-tangent dielectrics (Cirlex, Si, Ge)
are implemented but not yet placed in any geometry.

As of 2026-08-12 the `fix/core-hardening` and `light-pipe-example` branches are
**merged into `main`**, so everything below ships from a single branch: the core
simulation, the navigator/cache hardening, and the `BBRLightPipe` executable
with CAD (`.STL`) import.

### Recent correctness hardening

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
  x-face and sent the photon through solid copper (`crack_wall.mac`).
- `abspoints.term_status` is now `WorldExit` for photons killed at the world
  boundary; it used to echo the boundary process's stale per-thread status,
  including `BBRAbsorb` for photons that were never absorbed (`world_exit.mac`).
- The Planck emitter box is configurable (`/bbr/thermal/emitterCenter`,
  `/bbr/thermal/emitterSize`); `lightpipe.mac` sizes it to the bore, so photons
  are no longer created inside the light-pipe wall.
- Raw HFSS power ratios above 1 are capped at load time and exit-face-plane
  far-field directions get zero sampling weight; crack transmittance at normal
  incidence is 50.0% / 50.3% against the 50% ideal (`crack_transmit.mac`).

These changes were smoke-tested with a 10,000-event fixed-gun run using 15
workers; no geometry-navigation warnings, boundary-process errors, or stuck
tracks were observed. There is no registered `ctest` suite — correctness is
checked by the `scripts/check_*.py` PASS/FAIL validators, and
`scripts/run_regression.sh [build-dir]` builds, runs the eight regression macros
and executes every validator in one command (exit code = unexpected
failures; known reds are listed as XFAIL inside the script).

### Working

- **HFSS diffraction** — `BBRHFSSData` loads far-field + waveguide CSVs for one
  frequency; `BBRCrackLibrary` discovers each crack's frequency grid from the
  dataset directory names and lazy-loads the grid point nearest the photon; `BBSimOpBoundaryProcess` intercepts photons entering
  `vacuum_wg` crack volumes and routes them through the HFSS lookup.
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
  per-surface emissivity weighting. Default mode of `BBRTestPGA`; temperature
  set at runtime via `/bbr/thermal/setT`.

- **Test geometry** — `BBRTestDetectorConstruction`: 50 cm world, 4 mm Cu slab,
  two `vacuum_wg` crack daughters at z = 0 and z = 3 mm. Cu material is
  configurable via `/bbr/det/` before `/run/initialize`; gun and emitter
  settings via `/bbr/gun/` and `/bbr/thermal/` at any time.

- **Loss-tangent dielectrics** — `BBRMaterials::GetCirlex()`,
  `GetSiliconCrystal()`, `GetGermaniumCrystal()`: flat `RINDEX` plus
  `ABSLENGTH = c/(2πν·n·tanδ)` on the same 10 GHz–20 THz grid, with bulk
  absorption via stock `G4OpAbsorption` and Fresnel via the stock boundary
  process. Compile- and link-verified only — not yet placed in any geometry
  and not covered by a `check_*` validator.

- **ROOT output + analysis layer** — `BBRRunAction` / `BBRTestSteppingAction`
  write `output/bbr.root` (`crossings` + `abspoints` ntuples) and
  `output/bbr_legend.json` via `G4AnalysisManager` with ntuple merging under
  multithreading. `notebooks/copper_reflectance.ipynb` walks the copper Drude
  model end to end against the same module. `analysis/bbrsim/` (`io.py`,
  `physics.py`, `select.py`, `hfss.py`) is
  the single source of truth for loading and for the Drude / Planck /
  Hagen-Rubens formulas; every `check_*` and `plot_*` script reads through it.

### Not yet implemented

- Patched `G4OpBoundaryProcess` with `REFLECTIVITY` on `G4MaterialPropertiesTable`
  (upstream Geant4 PR target — BBRsim currently *wraps* the stock process rather
  than patching it)
- Leakage-current post-processing (the `abspoints` ntuple that feeds it exists;
  the loss-tangent-weighting analysis is Phase B2)
- PCB material; geometry placement and a transmission validator for the
  existing Si/Ge/Cirlex getters
- BBR calibration geometry (BB source + mesh-TES detector models)
- Anomalous-skin-effect correction to the classical Drude boundary optics

### Known scope limits

- **HFSS data is single-frequency; the lookup is not.** The dataset is chosen
  per photon by nearest frequency in log space, and the choice is recorded in
  the `hfss_freq_GHz` output column. Only a 500 GHz dataset exists per crack,
  so in practice every photon still gets the 500 GHz tables and broadband crack
  results stay indicative rather than quantitative. Dropping in real exports
  means adding `<id>_<freq>GHz_Ephi=N` directories; there is no interpolation
  between grid points. The mechanism is validated against a mock five-frequency
  tree (`crack_frequency.mac`).
- **Oblique incidence is validated at 45° only.** `crack_oblique.mac` and
  `check_crack_oblique.py` verify the azimuth fold/unfold on the θ = 135° row
  of the HFSS grid (the only oblique incidence in the data): momentum along
  the plate's long axis keeps its sign, the polarization filter behaves (E
  across the gap transmits, in-plane E is cut off), ±y and ±z tilts mirror.
  The z-mirror is unobservable on a parallel plate, and no HFSS run outside
  the [0°, 90°] wedge exists yet.

## Build

**Prerequisites:** Geant4 11.x (built with optical physics), CMake ≥ 3.16,
a C++17 compiler. Python scripts require the `bbrsim` conda environment
(NumPy, SciPy, Matplotlib, Pandas).

```bash
cmake --preset clang-release          # Apple Clang, Release, binaryDir build/
cmake --build --preset clang-release
cmake --install build                 # library + HFSS data into install/
```

The install step is required: with `BBRSIMDATA` unset and no `/bbr/dataDir`,
the executables read the HFSS data from `<install prefix>/share/BBRsim/data`
(the preset's prefix is `install/`), so a crack macro run before installing
aborts with `BBR011`.

`CMakePresets.json` pins `/usr/bin/clang`/`clang++` and `Release`; IDE CMake
integrations (VSCode CMake Tools) pick the preset up automatically instead of
offering their own compiler kits. The manual equivalent is
`cmake -S . -B build -DCMAKE_C_COMPILER=/usr/bin/clang
-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release`; a
`clang-debug` preset builds into `build-debug/`.

On macOS, configure with Apple Clang explicitly as shown. A bare `cmake ..` can
pick up Homebrew GCC (libstdc++), which compiles but fails at link against a
Geant4 built with Apple Clang (libc++) on every API whose signature contains
`std::` types. If you see undefined Geant4 symbols at link time, check
`CMAKE_CXX_COMPILER` in `build/CMakeCache.txt`.

Batch-only (no UI/visualization):

```bash
cmake -DWITH_GEANT4_UIVIS=OFF ..
make
make install
```

CMake copies all `.mac` files to the build directory alongside the executable.

## Run

The executable is `BBRSim` (not `OpNovice2`). All commands below run from the
`build/` directory.

### Planck emitter + crack diffraction

```bash
./BBRSim planck.mac          # 10 000 thermal photons at 4 K
./BBRSim test.mac            # 5M-event production run at 4 K
./BBRSim test_10K.mac        # 1M events at 10 K
```

Photons entering the `vacuum_wg` cracks are routed through the HFSS lookup;
expected transmittance at 500 GHz normal incidence is 50% per crack (the
unpolarized TEM-only ideal; observed 50.0% / 50.3%, see `crack_transmit.mac`
and `check_crack_transmittance.py`). Output is a ROOT file
`output/bbr.root` — two ntuples, `crossings` (one row per optical-photon
boundary crossing) and `abspoints` (one row per photon termination) — plus a
`output/bbr_legend.json` sidecar that maps the integer code columns
(status / event_type / volume / material) back to names. Runs are multithreaded;
G4Analysis merges the per-thread ntuples. Read it in Python via the shared
loader `analysis/bbrsim/io.py`, which decodes the codes to the legacy column
names:

```bash
conda run -n bbrsim python scripts/check_planck_spectrum.py build/output/bbr.root --temp 4
```

Every `check_*` / `plot_*` script reads BBRsim output through `analysis/bbrsim`;
none of them read a BBRsim-produced CSV. The only remaining `pandas.read_csv`
calls load external reference data (HFSS far-field tables and the Palik / Serov /
Geant4-IR copper comparison sets).

### Cu reflectance smoke test

```bash
./BBRSim reflectance.mac     # OFHC_Cu (RRR=100), 500 GHz, 10 000 events
./BBRSim test_of_cu.mac      # OF_Cu   (RRR=3),   500 GHz, 2 000 events
./BBRSim test_hp_cu.mac      # HP_Cu   (RRR=6),   500 GHz, 2 000 events
```

Output: `[BBR] reflectance mat=... N=... A_obs=... R_theory=...` printed to stdout.
At 4 K OFHC Cu sits on the relaxation plateau (D ≈ 4.9×10⁻⁵ at 500 GHz).
Compare against Drude theory:

```bash
./BBRSim reflectance.mac
conda run -n bbrsim python scripts/check_reflectance.py
```

Plot reflectance vs frequency across Cu grades and temperatures:

```bash
conda run -n bbrsim python scripts/plot_cu_reflectance.py
# output: build/cu_reflectance_plots.png
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

### Wrapper regression

`BBSimOpBoundaryProcess` falls through to the stock `G4OpBoundaryProcess` for
any geometry without `vacuum_wg` or `REFLECTIVITY` materials. The historical
fixed-seed regression macro (`verify_wrapper.mac`) was retired with the
OpNovice2 geometry split, and the stock OpNovice2 sources now live (unbuilt)
under `reference/opnovice2/`. Re-establishing an automated pass-through
regression — using a purpose-built minimal geometry that exercises stock
boundary optics, not the unbuilt OpNovice2 example — is open work. Until then,
when touching the wrapper, run the smoke tests (`reflectance.mac`,
`planck.mac`) and their `check_*` scripts before merging.

Five fixed-gun regression macros cover the wrapper's own physics:
`crack_wall.mac` (a photon inside a crack striking its Cu wall must reflect
about the wall normal and never enter the copper), `world_exit.mac` (photons
leaving the world are labelled `WorldExit` in `abspoints`), and
`crack_transmit.mac` (40k photons at normal incidence into crack1 must transmit
50% and never leave tangentially). Validate with
`check_crack_wall_reflection.py`, `check_no_photons_in_metal.py`,
`check_term_status.py` and `check_crack_transmittance.py`; the first two
invariant checks are meaningful on any output. `crack_oblique.mac` (16 runs at
45° off normal, random and fixed polarization via `/bbr/gun/pol`, one output
file per run) is validated by `check_crack_oblique.py` against the Python
mirror of the HFSS sampler in `analysis/bbrsim/hfss.py`. `crack_frequency.mac`
(16 runs spanning and straddling a mock five-frequency grid) is validated by
`check_crack_frequency.py`; it needs the mock tree first:

```bash
conda run -n bbrsim python scripts/make_mock_hfss_frequencies.py \
    --src data/waveguides --dst mock_hfss \
    --ids InfParallelPlate_crack1Rohan InfParallelPlate_crack2
cd build && ./BBRSim crack_frequency.mac
conda run -n bbrsim python scripts/check_crack_frequency.py build/output --data-dir mock_hfss
```

### Interactive (UI + visualization)

```bash
./BBRSim
```

Requires a Geant4 build with UI and visualization drivers enabled.

## Light-pipe geometry

A second executable, `BBRLightPipe`, models a 4 K → mixing-chamber light pipe.
It shares the physics list, materials, emitter, and ROOT output with `BBRSim`;
only the detector construction differs. Both executables are built by the
standard build above.

```bash
cd build && ./BBRLightPipe lightpipe.mac
```

`lightpipe.mac` places the shared Planck emitter just upstream of the warm
aperture and inside the 5 mm bore (`/bbr/thermal/emitterCenter -51 0 0 mm`,
`/bbr/thermal/emitterSize 1 7 7 mm`). With the `BBRSim` default emitter (a
1×20×20 mm patch centred at x = −50 mm) the emitting face sits inside the tube
wall and photons are created in the copper; `BuildParametric` warns (`LP002`)
if the configured emitter reaches into the wall.

Two build modes, selected by `/bbr/lightpipe/mode`:

- **`parametric`** (default) — a tube generated from bore radius, length, and
  wall thickness.
- **`cad`** — an ASCII `.STL` imported through the bundled header-only CADMesh
  (`library/include/CADMesh.hh`). The built-in reader is ASCII-only; binary STL needs
  assimp. A sample mesh ships at `data/cad/box_sample.stl`.

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

## License

Portions derived from the Geant4 `examples/extended/optical/OpNovice2` example
retain Geant4's license terms. BBRsim-specific additions will be released under
a compatible open-source license with the final Geant4 module publication.
