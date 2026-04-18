# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
mkdir build && cd build
cmake ..
make
```

To build without UI/visualization (batch-only):
```bash
cmake -DWITH_GEANT4_UIVIS=OFF ..
make
```

The CMake build automatically copies all `.mac` macro files to the build directory — the executable requires them to be co-located.

## Running the Simulation

**Interactive mode** (launches Geant4 UI/visualization):
```bash
cd build
./OpNovice2
```

**Batch mode** (run a specific macro):
```bash
cd build
./OpNovice2 electron.mac
./OpNovice2 boundary.mac
./OpNovice2 fresnel.mac
./OpNovice2 wls.mac
```

There is no automated test suite. Correctness is verified by comparing console output and histogram ROOT files against the reference `.out` files in the repository root.

## Architecture Overview

This is a Geant4 optical photon transport simulation studying Cerenkov radiation, scintillation, wavelength shifting (WLS/WLS2), and optical boundary interactions. The executable is named `OpNovice2`.

**Geometry:** World box (10m³, air) containing a 1m³ cubic tank (default: water). The configurable optical surface between the tank and world is the primary object of study.

**Physics stack:** FTFP_BERT base + G4EmStandardPhysics_option4 (replaces default EM) + G4OpticalPhysics. All configured in `OpNovice2.cc`.

**Data flow:**
- `PrimaryGeneratorAction` shoots particles (default: 3 eV optical photon at origin)
- `SteppingAction` records per-step surface interactions and optical process outcomes
- `Run` (custom G4Run subclass) accumulates statistics: Cerenkov/scintillation/WLS counts and energies, boundary process outcomes
- `HistoManager` books 26 predefined 1D histograms and writes ROOT/HDF5/XML/CSV output (default file: `opnovice2`)

**Interactive command system:** Three messenger classes expose runtime-configurable parameters:
- `DetectorMessenger` — surface finish, type, model, sigma_alpha, material properties (RINDEX, ABSLENGTH, etc.)
- `PrimaryGeneratorMessenger` — particle type, energy, position, direction, polarization
- `SteppingMessenger` — verbosity and stepping behavior

**Histogram layout:** Histos 1–9 cover optical spectra and timing (Cerenkov, scintillation, WLS); histos 10–26 cover boundary process outcomes and momentum directions per surface type.

**Key macro files:**
| Macro | Purpose |
|---|---|
| `electron.mac` | Cerenkov + scintillation from electron beam |
| `boundary.mac` | Sweeps all surface types and models |
| `fresnel.mac` | Fresnel reflection/transmission vs. angle |
| `wls.mac` | WLS/WLS2 processes |
| `coated.mac` | Thin-film coating effects |
| `complexRindex.mac` | Dielectric-metal surfaces with complex refractive index |
| `scint_by_particle.mac` | Particle-specific scintillation yields |
