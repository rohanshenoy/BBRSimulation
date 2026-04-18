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

## Current repository state (2026-04-18)

The repo is a thin BBR-oriented scaffold **on top of** the stock Geant4
`optical/OpNovice2` example. This is honest — most of the physics in the
list above is not yet implemented.

### BBR-specific additions

- `BBRSim.cc` — new entry point. Registers `BBSimPhysics` after
  `G4OpticalPhysics` so that `G4OpBoundaryProcess` already exists in the
  optical-photon process manager.
- `include/BBSimPhysics.hh`, `src/BBSimPhysics.cc` — `G4VPhysicsConstructor`
  subclass. `ConstructProcess()` locates `G4OpBoundaryProcess` in the
  optical-photon process manager, removes it, wraps it in
  `BBSimOpBoundaryProcess`, and re-registers the wrapper. Pattern follows
  SuperSim's `CDMSRDecayPhysics::WrapRDMProcess()`.
- `include/BBSimOpBoundaryProcess.hh`, `src/BBSimOpBoundaryProcess.cc` —
  `G4WrapperProcess` subclass. **Currently a pure pass-through** whose
  `PostStepDoIt` forwards verbatim to the wrapped `G4OpBoundaryProcess`.
  This is the injection point where custom BBR boundary physics
  (diffraction / vacuum→copper tabulated reflectance / HFSS interface)
  will live.
- `verify.mac` — fixed-seed regression macro. Output must be bit-identical
  to an unwrapped run until the wrapper stops being a pass-through.
- `CMakeLists.txt` — builds a single executable `BBRSim` from `BBRSim.cc`
  plus everything in `src/`.

### Stock OpNovice2 pieces (still unchanged)

- `OpNovice2.cc` — original entry point, retained for reference.
- `DetectorConstruction`, `PrimaryGeneratorAction`, `SteppingAction`,
  `Run`, `HistoManager`, `TrackingAction`, and all three messenger
  classes — unmodified OpNovice2. They shoot a single particle at a
  configurable 1 m³ box inside a 10 m³ world and histogram Cerenkov /
  scintillation / WLS / boundary-process outcomes.
- `.mac` files (`electron.mac`, `boundary.mac`, `fresnel.mac`, `wls.mac`,
  `coated.mac`, `complexRindex.mac`, `scint_by_particle.mac`,
  `OpNovice2.mac`, `vis.mac`) — all stock OpNovice2 test cases. **None
  exercise BBR physics.**

### What is explicitly *not* in the repo yet

- No `ThermalSurface` / `GetBBSpecCDF` / `BBEvt` / `GeometricSurface` —
  no way to emit a Planck-distributed photon from a surface today.
- No patched `G4OpBoundaryProcess` — `REFLECTIVITY` still lives only on
  `G4OpticalSurface`.
- No CADMesh integration, no `.STL` import.
- No HFSS `.csv` loader, no bounded-volume handoff in `SteppingAction`.
- No cryogenic-material optical-property database.
- No BBR calibration geometry (BB source or mesh-TES detector model).

Treat the repo as **scaffolding**: `BBSimOpBoundaryProcess` is the
verified hook point, but the physics that will eventually live behind it
is still to be written.

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
./BBRSim electron.mac        # stock OpNovice2 test case
./BBRSim verify.mac          # wrapper pass-through regression
```

The executable is `BBRSim`, **not** `OpNovice2`.

## Verification discipline

`BBSimOpBoundaryProcess` is a load-bearing hook: every subsequent piece of
BBR physics will be added behind it. While it is a pass-through, the
guarantee is **bit-identical output** against an unwrapped `G4OpBoundaryProcess`
run. `verify.mac` is the regression harness — run with fixed seeds and
compare histograms before merging anything that touches the wrapper or
`BBSimPhysics::WrapOpBoundaryProcess`.

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
