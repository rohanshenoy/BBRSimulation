# CLAUDE.md

Guidance for Claude / Claude Code when working on this repository.

## Project: BBRsim

This repository is the early-stage implementation of **BBRsim**, a Geant4-based
free-space blackbody radiation (BBR) simulation for cryogenic experiments. The
work continues Yen-Yung Chang's Caltech thesis (2023, Ch. 5) and is the subject
of the NSF QIS proposal led by Sunil Golwala (Caltech) / Nader Mirabolfathi
(TAMU), with Rohan Shenoy as the student lead.

### Scientific motivation

Sub-Kelvin devices — superconducting-circuit qubits (transmons, fluxonium) and
phonon-mediated dark matter / CEvNS detectors (SuperCDMS, TESSERACT) — are
limited by non-thermal populations of broken Cooper pairs / excited charge
carriers. A known contributor is free-space BBR leaking through machining gaps,
cable feedthroughs, and flange-lid joints in otherwise "sealed" cryostat
chambers. BBR from 4 K surfaces peaks near 1.3 mm; mating tolerances of
O(100 µm) leave every chamber a leaky waveguide for long-wavelength photons
(cutoff λ = 2 d_max; bends do **not** reflect waveguide modes). Ray-trace and
commercial non-sequential optics tools are both inadequate: BBR is emitted by
*every* warm surface and must be propagated as rays in open volumes **and** as
waves in narrow gaps. BBRsim's job is to do both, consistently, on CAD-accurate
geometry.

### Project objectives (from the NSF proposal)

- **O1** — Finish the free-space BBR simulation package (this repo).
- **O2** — Validate experimentally: a temperature-controlled TK-RAM / OFHC-Cu
  BB source (4–20 K) + purpose-built mesh-TES photon detectors, first at
  TAMU's DR, then on a SuperCDMS SNOLAB Pathfinder tower at SLAC.
- **O3** — Release BBRsim publicly as a Geant4 module; upstream the required
  `G4OpBoundaryProcess` modifications via the Geant4 EM Physics Working Group.

## Physics architecture (what BBRsim must do)

Two coupled propagation regimes, switched per-volume by the size criterion
in Chang Table 5.1 (roughly: a space is "open" if its smallest dimension is
≥10× the photon wavelength, else "bounded"):

1. **Open volumes → ray-trace (Geant4).** Built-in `G4OpticalPhoton` with
   all other processes disabled. Straight-line propagation with a
   frequency-dependent mean free path in dielectrics (Si, Ge, polymers) and
   customized reflectance at metal surfaces. Fresnel's equations are
   *unreliable* for cryogenic metals in the trans-mm band (reflectance is
   dominated by the imaginary refractive index, which varies orders of
   magnitude with RRR and surface prep), so BBRsim overrides boundary
   optics with tabulated per-material `REFLECTIVITY`.
2. **Bounded volumes → wave (pre-computed HFSS S-parameters).** For gaps
   with a dimension ≤ ~10× λ, the photon's reflection/transmission
   probabilities and outgoing k-direction PDFs come from ANSYS HFSS sweeps
   over (k_in, ν, polarization), stored as `.csv` and sampled at the
   interface. The bend geometry of a flange-lid joint or cable slot is
   solved once in HFSS; BBRsim then looks up the transmittance.

### Required new Geant4 components

- **Emission** — `ThermalSurface`, `GeometricSurface`, `GetBBSpecCDF`,
  `BBEvt` classes; Planck PDF sampler; per-surface emissivity; T³-weighted
  surface selection in a custom `PrimaryGeneratorAction`. Photons are
  emitted isotropically into the outward hemisphere with no polarization
  preference (valid after many reflections even if not at emission).
- **Modified `G4OpBoundaryProcess`** — `REFLECTIVITY` lifted from
  `G4OpticalSurface` onto `G4MaterialPropertiesTable`, so per-material
  optical properties can be set at CAD-import time without instantiating a
  logical surface per object. (Chang's thesis published a patched
  `G4OpBoundaryProcess.cc`; upstream release is an O3 deliverable.)
- **CAD import** — CADMesh (2nd ver.) imports `.STL` as `G4TessellatedMesh`
  and hierarchical `G4AssemblyVolume`. CAD comes from SOLIDWORKS at max
  resolution (~7 µm, 0.5°). Parts grouped into "sub-detector" C++ classes
  (cryostat, BB shield, detector module) for composability.
- **Open↔bounded interface** — custom hooks in `UserSteppingAction` that
  sample HFSS-tabulated reflection/transmission when a photon enters a
  flagged gap volume.
- **Material optical-property database** — one shared class; current
  materials targeted are vacuum, OFHC Cu, Cirlex (polymer proxy), PCB,
  crystalline Si/Ge, perfect reflector/absorber. Literature data on 20+
  other mm-wave absorbers is compiled for future use.
- **Absorption / leakage-current post-processing** — photon termination
  events in Si/Ge are recorded via a ROOT `TTree` (`ABSPoint` struct:
  energy, momentum, position, reflection count, terminating volume), and
  leakage-current rates are computed offline by weighting each absorption
  by the ratio of the shallow-impurity-only mean free path to the
  empirical loss-tangent-derived mean free path. Use `tan δ_Si = 1e-4`,
  `tan δ_Ge = 6e-5` as defaults (Chang §5.3.1.4).

## Current repository state (2026-05-05)

The repo has progressed beyond the initial scaffold. The HFSS diffraction
path is fully implemented and validated. A material framework (`BBRMaterials`,
`vacuum_wg`) is in place. The next two pieces — material reflectance tables
and the Planck emitter — are designed but not yet built.

### BBR-specific additions (implemented)

- `BBRSim.cc` — entry point. Routes to `BBRCrackDetectorConstruction` +
  `BBRDiffractionActionInit` when the mac name contains `diffraction` or
  `validation`; falls back to OpNovice2 geometry otherwise.
- `include/BBSimPhysics.hh`, `src/BBSimPhysics.cc` — `G4VPhysicsConstructor`
  that wraps `G4OpBoundaryProcess` with `BBSimOpBoundaryProcess`. Pattern
  follows SuperSim's `CDMSRDecayPhysics::WrapRDMProcess()`.
- `include/BBSimOpBoundaryProcess.hh`, `src/BBSimOpBoundaryProcess.cc` —
  `G4WrapperProcess` subclass. Intercepts photons that enter a `vacuum_wg`
  volume and routes them through `HandleDiffractionBoundary` (HFSS lookup).
  Everything else falls through to the stock `G4OpBoundaryProcess`.
- `include/BBRMaterials.hh` — static factory for `vacuum_wg` (a near-vacuum
  material with RINDEX=1 used to flag crack volumes). Seed of the planned
  `BBRMaterialDB`.
- `include/BBRHFSSData.hh`, `src/BBRHFSSData.cc` — loads `far_field.csv` +
  `waveguide.csv` for one HFSS dataset (indexed by `(IWavePhi, IWaveTheta)`
  angle grid). Provides `GetTransmittance`, `SampleOutgoingDirection`,
  `SampleExitPosition`. Runtime CDFs handle arbitrary polarization including
  cross terms.
- `include/BBRCrackLibrary.hh`, `src/BBRCrackLibrary.cc` — Meyer's singleton
  that lazy-loads `BBRHFSSData` per dataset ID (= `vacuum_wg` volume name,
  strip any `:N` suffix). Adding a new crack requires only placing a new
  `vacuum_wg` volume — no code changes.
- `include/BBRCrackDetectorConstruction.hh`, `src/BBRCrackDetectorConstruction.cc` —
  world + two crack volumes: `InfParallelPlate_crack1_500GHz` (b=50 µm, z=0)
  and `InfParallelPlate_crack2_500GHz` (b=100 µm, z=3 mm), both `vacuum_wg`.
- `include/BBRDiffractionPGA.hh`, `src/BBRDiffractionPGA.cc` — fires one
  optical photon per event at (−20 mm, 0, z) along +x; 500 GHz; polarization
  (0, −1/√2, −1/√2). Gun Z exposed via `/bbr/gun/setZ <value> mm` messenger.
- `include/BBRDiffractionActionInit.hh`, `src/BBRDiffractionActionInit.cc` —
  wires `BBRDiffractionPGA` + stepping action; accepts `gunZ_mm` at construction.
- `include/BBRDiffractionSteppingAction.hh`, `src/BBRDiffractionSteppingAction.cc` —
  records transmitted photon position + crack ID to `diffraction_output.csv`.
- `verify_wrapper.mac` — fixed-seed OpNovice2 regression (renamed from
  `verify.mac`). Bit-identical output required before/after any wrapper change.
- `diffraction.mac` — crack1 smoke test. T_obs ≈ 0.528, theory ≈ 0.527.
- `diffraction_crack2.mac` — crack2 smoke test, gun at z=3 mm.
- `validation.mac` — combined crack1 + crack2 in one session;
  `/bbr/gun/setZ 3 mm` between runs. Output in `diffraction_output.csv`
  with `crack_id` column. T_exp: crack1 ≈ 52.7%, crack2 ≈ 50.4%.
- `plot_diffraction.py` — 3-panel plot: per-crack transmittance vs. event
  count, exit-position distribution. Run with `conda run -n bbrsim python`.

### Stock OpNovice2 pieces (unchanged)

- `OpNovice2.cc` — original entry point, retained for reference.
- `DetectorConstruction`, `PrimaryGeneratorAction`, `SteppingAction`,
  `Run`, `HistoManager`, `TrackingAction`, and all three messenger
  classes — unmodified OpNovice2.
- `.mac` files (`electron.mac`, `boundary.mac`, `fresnel.mac`, `wls.mac`,
  `coated.mac`, `complexRindex.mac`, `scint_by_particle.mac`,
  `OpNovice2.mac`, `vis.mac`) — stock OpNovice2 test cases. **None
  exercise BBR physics.**

### What is explicitly *not* in the repo yet

- **No Planck emitter** — `BBRThermalPGA` / `BBRPlanckSampler` / `BBRThermalSurface`
  not yet built. All current runs fire a monochromatic 500 GHz test photon.
- **No material reflectance tables** — vacuum→OFHC Cu (Hagen-Rubens),
  perfect absorber, perfect reflector not yet implemented. `BBR_REFLECTIVITY`
  on `Material2`'s MPT is the planned injection point in `BBSimOpBoundaryProcess`.
- No patched `G4OpBoundaryProcess` — `REFLECTIVITY` still lives only on
  `G4OpticalSurface` in stock Geant4 11.4.
- No CADMesh integration, no `.STL` import.
- No cryogenic-material optical-property database beyond `vacuum_wg`.
- No ROOT `TTree` output or leakage-current post-processing.
- No BBR calibration geometry (BB source or mesh-TES detector model).

## Implementation log

Plans are stored in `~/.claude/plans/`. Completed plans are listed here for
reference; do not re-implement work that is already in the git history.

### [DONE] Pass-through wrapper skeleton
`majestic-floating-reef.md` — `BBSimPhysics` + `BBSimOpBoundaryProcess`
(pure pass-through). Verified bit-identical to unwrapped `G4OpBoundaryProcess`
via `verify_wrapper.mac`.

### [DONE] HFSS diffraction physics
`scalable-wandering-shore.md` — `BBRHFSSData`, `HandleDiffractionBoundary`,
coordinate-frame convention (HFSS ↔ Geant4 world). Validated T_obs ≈ 0.528
at 500 GHz normal incidence (theory 52.7%).

### [DONE] Generalized crack detection (`vacuum_wg` material)
`dreamy-hopping-llama.md` — crack volumes identified by `vacuum_wg` material
(not magic string). Volume name = HFSS dataset ID. Orientation extracted from
touchable rotation at runtime.

### [DONE] BBRCrackLibrary + crack2 geometry + combined validation macro
`05072026-crack-library.md` (archived) — `BBRCrackLibrary` singleton owns
dataset cache. crack2 placed at z=3 mm. `validation.mac` runs both cracks
in one session via `/bbr/gun/setZ`. Per-crack CSV output + 3-panel plots.

## Upcoming work

Plans live in `docs/superpowers/plans/` (gitignored). Naming: `mmddyyyy-feature.md`.

### [BLOCKED — awaiting Yen-Yung Chang's code] Material reflectance + Planck emitter
`docs/superpowers/plans/05072026-reflectance-and-planck.md` —
Full implementation plan covering both features. Decision: specular-only
reflection (Rayleigh criterion: surface roughness << λ at trans-mm).

**Reflectance scope:** `OFHC_Cu` material with Hagen-Rubens `BBR_REFLECTIVITY`
table (20 log-spaced energies, 50 GHz–20 THz, σ=5.96×10⁹ S/m for RRR=100);
`BBR_PerfectAbsorber`; `BBR_PerfectReflector`. Intercept in
`BBSimOpBoundaryProcess::PostStepDoIt` before stock fall-through.
New: `BBRReflectanceDetectorConstruction`, `reflectance.mac`.

**Planck emitter scope:** `BBRPlanckSampler` (header-only CDF table for
u³/(e^u−1)), `BBRThermalSurface` (POD struct), `BBRThermalPGA` (hardcoded
4 K patch), `BBRPlanckActionInit`, `BBRPlanckDetectorConstruction`,
`planck.mac`, `scripts/check_planck_spectrum.py`.

## Build

```bash
mkdir build && cd build
cmake ..
make
```

Batch-only (no UI/vis):
```bash
cmake -DWITH_GEANT4_UIVIS=OFF ..
make
```

CMake copies all root-level `.mac` files next to the built executable.

## Run

Interactive (UI + visualization):
```bash
cd build
./BBRSim
```

Batch (explicit macro):
```bash
./BBRSim electron.mac          # stock OpNovice2 test case
./BBRSim verify_wrapper.mac    # wrapper pass-through regression
./BBRSim diffraction.mac       # crack1 HFSS diffraction smoke test
./BBRSim validation.mac        # crack1 + crack2 combined validation
./BBRSim reflectance.mac       # OFHC Cu reflectance smoke test (planned)
./BBRSim planck.mac            # Planck emitter test, writes planck_output.csv (planned)
```

The executable is `BBRSim`, **not** `OpNovice2`.

## Verification discipline

`BBSimOpBoundaryProcess` is a load-bearing hook: every subsequent piece of
BBR physics will be added behind it. `verify_wrapper.mac` is the regression
harness — run with fixed seeds and compare histograms before merging anything
that touches the wrapper or `BBSimPhysics::WrapOpBoundaryProcess`. Geometries
without `vacuum_wg` volumes fall through identically to stock `G4OpBoundaryProcess`.

When new physics is added behind the wrapper (e.g. a tabulated vacuum→Cu
reflectance or an HFSS waveguide handoff), it must be selectable via a
messenger or surface flag so the pass-through path remains testable.

## References

- Chang, Y.-Y. (2023). *SuperCDMS HVeV Run 2 Low-Mass Dark Matter Search,
  Highly Multiplexed Phonon-mediated Particle Detector with Kinetic
  Inductance Detector, and the Blackbody Radiation in Cryogenic
  Experiments*. PhD thesis, Caltech. **Chapter 5** is the primary
  physics + simulation reference for this project.
- Golwala & Mirabolfathi, *BBRsim: Free-Space Blackbody Radiation
  Simulation for Superconducting Circuits and Cryogenic Detectors*, NSF
  QIS proposal (2025). Defines O1–O3 and the validation plan.
- Agostinelli et al., "GEANT4 — a simulation toolkit," NIM A **506**, 250
  (2003).
- Poole et al., CADMesh (2nd ver.) — CAD → Geant4 geometry import.

## Working notes for AI coding agents

- The executable is `BBRSim`. Do not reintroduce `./OpNovice2` in
  documentation, tooling, or CMake.
- The `.mac` files inherited from OpNovice2 are **regression inputs**, not
  representative BBR workloads. Do not treat their output as evidence of
  BBR physics working.
- There is no automated test suite. Correctness for the wrapper is
  verified by diffing fixed-seed histogram output (`verify.mac` →
  `verify.root`) against an unwrapped baseline.
- When adding new BBR physics, prefer new classes (e.g.
  `BBRThermalSurface`, `BBRMaterialDB`) over editing the OpNovice2
  `DetectorConstruction` / `PrimaryGeneratorAction` in place — the stock
  files are kept around as a known-good reference.
- The planned patch to `G4OpBoundaryProcess.cc` (moving `REFLECTIVITY`
  onto `G4MaterialPropertiesTable`) is upstream-facing work; any local
  implementation should be isolated so it can be turned into a Geant4
  PR cleanly.
- Python analysis scripts (`plot_diffraction.py`, future `plot_planck.py`)
  must be run as `conda run -n bbrsim python <script>`. Plain `python3`
  does not have the required packages even if they appear installed.
