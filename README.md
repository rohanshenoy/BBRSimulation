# BBRsim

Geant4-based free-space blackbody radiation (BBR) simulation for cryogenic
experiments. Continues the work of Yen-Yung Chang (Caltech PhD, 2023, Ch. 5)
and implements the technical scope of the Golwala / Mirabolfathi NSF QIS
proposal *BBRsim: Free-Space Blackbody Radiation Simulation for
Superconducting Circuits and Cryogenic Detectors* (2025).

## Why this exists

Superconducting-circuit qubits and sub-Kelvin phonon-mediated particle
detectors are limited by non-thermal quasiparticle / charge-carrier
populations. One cause is free-space BBR that leaks through machining
gaps, cable slots, and flange-lid joints in otherwise "sealed" cryostat
chambers. BBR from a 4 K surface peaks near λ = 1.3 mm, comparable to or
larger than typical mating tolerances, so every nominally closed chamber
behaves as a leaky waveguide for long-wavelength photons. Existing tools
do not handle this well: ray-trace simulators assume open propagation,
and non-sequential optical simulators assume discrete sources rather than
an "everywhere" thermal source.

BBRsim's goal is to do both — ray-trace in open volumes, wave propagation
(via pre-computed HFSS S-parameters) in gaps — on CAD-accurate geometry,
so apparatus can be designed *against* BBR backgrounds instead of
debugging them afterward.

## Status

**Early scaffolding on top of Geant4's `optical/OpNovice2` example.**

What exists:

- `BBRSim.cc` — executable entry point.
- `BBSimPhysics` — `G4VPhysicsConstructor` that wraps the
  `G4OpBoundaryProcess` registered by `G4OpticalPhysics`.
- `BBSimOpBoundaryProcess` — `G4WrapperProcess` subclass; pure
  pass-through today. This is the hook point where BBR-specific boundary
  physics (tabulated cryogenic metal reflectance, HFSS waveguide
  transmission, diffraction) will be injected.
- `verify.mac` — fixed-seed regression to prove the wrapper stays
  bit-identical to an unwrapped run while it is a pass-through.
- Unmodified OpNovice2 `DetectorConstruction`, `PrimaryGenerator`,
  `SteppingAction`, `Run`, `HistoManager`, messengers, and `.mac` test
  cases — retained as a known-good baseline.

What is **not** here yet (all planned): a Planck-distributed
`ThermalSurface` emission system; the patched `G4OpBoundaryProcess` that
moves `REFLECTIVITY` from `G4OpticalSurface` to
`G4MaterialPropertiesTable`; CADMesh / SOLIDWORKS `.STL` import; an HFSS
`.csv` loader and open↔bounded interface in `UserSteppingAction`; a
cryogenic-material optical-property database; and the BB source +
mesh-TES detector geometries needed for validation.

`CLAUDE.md` has a more detailed breakdown of the physics architecture
and the implementation gap.

## Build

```bash
mkdir build && cd build
cmake ..
make
```

Batch-only (no UI/visualization):

```bash
cmake -DWITH_GEANT4_UIVIS=OFF ..
make
```

CMake copies all `.mac` files to the build directory alongside the
executable.

## Run

Interactive (UI + vis):

```bash
cd build
./BBRSim
```

Batch:

```bash
./BBRSim electron.mac     # stock OpNovice2 — Cerenkov + scintillation from e-
./BBRSim boundary.mac     # stock OpNovice2 — boundary process sweep
./BBRSim verify.mac       # fixed-seed wrapper regression
```

The executable is `BBRSim`. (The original OpNovice2 `main()` is also
present as `OpNovice2.cc` for reference; it is not built by the default
`CMakeLists.txt`.)

## Project plan and team

From the NSF proposal, period of performance 2026Q4–2029Q3:

- **O1** — Complete BBRsim development.
- **O2** — Validate against a temperature-controlled TK-RAM / OFHC-Cu BB
  source (4–20 K, ~160 cm² emitting area, CFRP thermal break) plus
  purpose-built mesh-TES photon detectors (W TES, 3 mm on Si, three
  designs optimized for 4 / 9 / 20 K peak emission). First stages use
  manufactured test geometries at TAMU; the final stage uses a
  SuperCDMS SNOLAB Pathfinder detector tower at SLAC.
- **O3** — Publish and release BBRsim as a Geant4 module; upstream the
  `G4OpBoundaryProcess` `REFLECTIVITY` change via the Geant4 EM Physics
  Working Group.

Team: Golwala (PI, Caltech), Mirabolfathi (co-PI, TAMU), Shenoy (student
lead, Caltech), Xiong (Brinson postdoc, Caltech), Kurinsky / Partridge
(SLAC SuperCDMS, unfunded collaborators), Chang (BBRsim originator,
unfunded advisor).

## References

- Chang, Y.-Y., *SuperCDMS HVeV Run 2 Low-Mass Dark Matter Search,
  Highly Multiplexed Phonon-mediated Particle Detector with Kinetic
  Inductance Detector, and the Blackbody Radiation in Cryogenic
  Experiments* (Caltech PhD thesis, 2023). Chapter 5 is the primary
  reference for BBRsim physics and the particle-like simulation design.
- Golwala & Mirabolfathi, *BBRsim* NSF QIS proposal (2025).
- Agostinelli et al., "GEANT4 — a simulation toolkit," *NIM A* **506**,
  250 (2003).
- Poole et al., CADMesh (2nd ver.) — CAD → Geant4 geometry import.

## License

Portions derived from the Geant4 `examples/extended/optical/OpNovice2`
example retain Geant4's license terms. BBRsim-specific additions will be
released under a compatible open-source license with the final Geant4
module publication.
