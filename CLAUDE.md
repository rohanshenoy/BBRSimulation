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
- **Material optical-property database** — one shared class (`BBRMaterials`
  namespace, header-only). Metals use the full **Drude model** (Griffiths
  §9.4) parameterized by (RRR, T_K); Hagen-Rubens is only the low-frequency
  limit (ωτ << 1) and breaks down above ~64 GHz for OFHC Cu at 4 K. For
  copper: `σ_DC(T) = RRR × σ_RT` where σ_RT = 5.96×10⁷ S/m is universal;
  τ = σ_DC × mₑ/(nₑ e²); R = |(ñ−1)/(ñ+1)|² from the full complex
  refractive index. Current materials: vacuum, OFHC Cu, Cirlex (polymer
  proxy), PCB, crystalline Si/Ge, perfect reflector/absorber. See
  `docs/physics/copper_reflectance_model.md` for the full derivation.
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

### [UNBLOCKED] Drude copper reflectance model

Replace `BuildHagRubMaterial` (Hagen-Rubens) with `BuildDrudeMaterial` in
`include/BBRMaterials.hh`. Parameters: `(name, RRR, T_K)`. Derives σ_DC via
Matthiessen's rule, τ from the Drude formula, then evaluates full Griffiths
§9.4 reflectance at each frequency. Replaces all three named Cu getters with
a single `GetCopper(RRR, T_K)`. Update `GetCopperByName()` messenger to
accept RRR as integer and optional temperature. See
`docs/physics/copper_reflectance_model.md`.

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

### Copper reflectance physics

- σ_RT = 5.96×10⁷ S/m — universal for all Cu grades at 273 K; does **not**
  depend on material quality.
- σ(4 K) = RRR × σ_RT. RRR range: 1 (disordered) → ~10 (commercial) →
  ~100 (OFHC) → ~500 (ultra-pure crystal).
- τ = σ_DC × mₑ / (nₑ e²). Cu constants: nₑ = 8.49×10²⁸ m⁻³,
  mₑ = 9.109×10⁻³¹ kg.
- Hagen-Rubens valid only for f << 1/(2πτ). OFHC Cu at 4 K (RRR=100):
  τ ≈ 2.5 ps → H-R valid below ~64 GHz. BBRsim starts at 50 GHz — always
  use the full Drude model.
- Matthiessen's rule: 1/σ(T) = 1/σ_impurity + 1/σ_phonon(T), where
  1/σ_impurity = 1/(RRR × σ_RT) and σ_phonon(T) ≈ σ_RT × 273/T (T > 50 K).
- Full Drude: σ(ω) = σ_DC/(1−iωτ); insert into Griffiths §9.4 to get k̃;
  R = |(ñ−1)/(ñ+1)|² where ñ = c·k̃/ω.
- **Do NOT hardcode σ per Cu variant.** Derive everything from (RRR, T_K).
  See `docs/physics/copper_reflectance_model.md`.

## Geant4 API Reference

This section is a practical reference for the Geant4 11.x APIs used in BBRsim,
derived from the Geant4 Beginner Course and BBRsim's own codebase.

### Mandatory initialization pattern

Every application registers exactly three things with the run manager:

```cpp
auto* runManager = G4RunManagerFactory::CreateRunManager(); // auto-selects MT/sequential
runManager->SetUserInitialization(new YourDetectorConstruction());
runManager->SetUserInitialization(physicsListPtr);          // G4VModularPhysicsList
runManager->SetUserInitialization(new YourActionInitialization());
// Then: runManager->Initialize(); runManager->BeamOn(N);
```

Use `G4PhysListFactory` to get a reference physics list by name (e.g. `"FTFP_BERT_EMZ"`).
For BBRsim, `BBSimPhysics` is a `G4VPhysicsConstructor` added to a modular list.

### Geometry: solid → logical → physical

```cpp
// 1. Solid (shape, no material, no position)
G4Box* solid = new G4Box("name", halfX, halfY, halfZ);   // args are HALF-lengths

// 2. Logical volume (solid + material, no position)
G4LogicalVolume* lv = new G4LogicalVolume(solid, material, "name");

// 3. Physical volume (logical + position inside a mother)
new G4PVPlacement(
    nullptr,               // rotation (G4RotationMatrix* or nullptr)
    G4ThreeVector(x,y,z),  // translation
    lv,                    // logical volume to place
    "name",                // physical volume name
    motherLV,              // mother logical volume (nullptr for world)
    false,                 // pMany (unused, always false)
    0);                    // copy number
```

World volume: mother = `nullptr`. `Construct()` must return the world
`G4VPhysicalVolume*`. Other common solids: `G4Tubs`, `G4Sphere`, `G4Trd`,
`G4SubtractionSolid` (CSG boolean).

Rotations:
```cpp
G4RotationMatrix* rot = new G4RotationMatrix();
rot->rotateX(90.*deg);   // then pass to G4PVPlacement
```

### Materials

**NIST database** (preferred for standard materials):
```cpp
G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
// Names: G4_Si, G4_Ge, G4_Cu, G4_STAINLESS-STEEL, G4_AIR, G4_WATER, etc.
```

**Custom material by hand:**
```cpp
G4Material* mat = new G4Material("name", Z, A_g_per_mole*g/mole,
                                  density*g/cm3, kStateSolid, temperature*kelvin);
// Multi-element:
G4Element* el = new G4Element("name","symbol", Z, A*g/mole);
G4Material* mat = new G4Material("name", density*g/cm3, nComponents);
mat->AddElement(el, massFraction);
```

**Material Properties Table** (required for optical photons):
```cpp
G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
// Constant property:
mpt->AddConstProperty("RINDEX", 1.0);
// Energy-indexed vector (energies in ASCENDING order):
std::vector<G4double> energies = {1.0*eV, 2.0*eV, 3.0*eV};
std::vector<G4double> rindex   = {1.5,    1.5,    1.5};
mpt->AddProperty("RINDEX", energies, rindex);
mat->SetMaterialPropertiesTable(mpt);
```

Key optical property names on `G4Material` MPT: `RINDEX`, `ABSLENGTH`,
`SCINTILLATIONYIELD`, `RAYLEIGH`. BBRsim adds `BBR_REFLECTIVITY` as a
custom property read by `BBSimOpBoundaryProcess`.

### Optical surfaces and boundary process

```cpp
// Logical surface between two logical volumes:
G4OpticalSurface* opSurf = new G4OpticalSurface("name");
opSurf->SetType(dielectric_metal);   // dielectric_metal, dielectric_dielectric, etc.
opSurf->SetModel(unified);           // unified, glisur, LUT, davis
opSurf->SetFinish(polished);         // polished, ground, polishedfrontpainted, etc.

G4MaterialPropertiesTable* surfMPT = new G4MaterialPropertiesTable();
surfMPT->AddProperty("REFLECTIVITY", energies, reflectivities);
opSurf->SetMaterialPropertiesTable(surfMPT);

// Attach to a pair of volumes (border surface: directed vol1→vol2):
new G4LogicalBorderSurface("name", physVol1, physVol2, opSurf);
// OR skin surface (wraps an entire logical volume):
new G4LogicalSkinSurface("name", logVol, opSurf);
```

`G4OpBoundaryProcess` runs automatically when `G4OpticalPhoton` crosses a
boundary. It checks for a `G4LogicalBorderSurface` or `G4LogicalSkinSurface`
first; if none, it uses the `RINDEX` of both materials to compute Fresnel.

BBRsim overrides this via `BBSimOpBoundaryProcess` (`G4WrapperProcess`
subclass). The wrapper intercepts `PostStepDoIt`, checks if the volume
material is `vacuum_wg`, and routes to HFSS lookup before (or instead of)
calling the stock process.

### Physics list / constructor pattern

```cpp
// G4VPhysicsConstructor subclass:
class BBSimPhysics : public G4VPhysicsConstructor {
public:
    void ConstructParticle() override { /* define particles if needed */ }
    void ConstructProcess() override  {
        // Wrap or register processes:
        WrapOpBoundaryProcess();
    }
private:
    void WrapOpBoundaryProcess();
};

// Registering in a modular list:
auto* pl = new G4VModularPhysicsList();
pl->RegisterPhysics(new G4OpticsPhysics());   // stock optical
pl->RegisterPhysics(new BBSimPhysics());      // BBR-specific override
runManager->SetUserInitialization(pl);
```

### G4WrapperProcess pattern

```cpp
class BBSimOpBoundaryProcess : public G4WrapperProcess {
public:
    G4VParticleChange* PostStepDoIt(const G4Track& track,
                                     const G4Step& step) override {
        // intercept:
        if (IsVacuumWG(track)) return HandleDiffractionBoundary(track, step);
        // fall through to stock:
        return G4WrapperProcess::PostStepDoIt(track, step);
    }
};

// Wrapping the existing process (in BBSimPhysics::ConstructProcess):
G4ProcessManager* pm = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
G4int idx = pm->GetProcessIndex(stockBoundaryProc, idxPostStep);
auto* wrapper = new BBSimOpBoundaryProcess();
wrapper->RegisterProcess(stockBoundaryProc);   // store original inside wrapper
pm->RemoveProcess(stockBoundaryProc);
pm->AddDiscreteProcess(wrapper);
```

### Action initialization (MT-safe pattern)

```cpp
class YourActionInitialization : public G4VUserActionInitialization {
public:
    // Called once by master thread only — register RunAction for master:
    void BuildForMaster() const override {
        SetUserAction(new YourRunAction());
    }
    // Called by each worker thread — register all actions:
    void Build() const override {
        auto* pg = new YourPrimaryGeneratorAction();
        SetUserAction(pg);
        auto* run = new YourRunAction(pg);
        SetUserAction(run);
        auto* evt = new YourEventAction();
        SetUserAction(evt);
        SetUserAction(new YourSteppingAction(evt));
    }
};
```

Rule: objects created in `Build()` are **worker-thread-local**. Never share
mutable state between them without a mutex. The detector construction pointer
is safe to pass down (it's read-only after `Construct()` returns).

### Stepping action

```cpp
void YourSteppingAction::UserSteppingAction(const G4Step* step) {
    // Volume check:
    auto* preVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    // Per-step data:
    G4double eDep   = step->GetTotalEnergyDeposit();
    G4double trackL = step->GetStepLength();

    // Track / particle info:
    const G4Track* track = step->GetTrack();
    const G4ParticleDefinition* pDef = track->GetParticleDefinition();
    G4ThreeVector pos  = track->GetPosition();
    G4ThreeVector dir  = track->GetMomentumDirection();
    G4double      eKin = track->GetKineticEnergy();

    // Kill a track:
    track->SetTrackStatus(fStopAndKill);  // via non-const pointer

    // Post-step point (after the step):
    auto* post = step->GetPostStepPoint();
    G4String procName = post->GetProcessDefinedStep()->GetProcessName();
}
```

### Run accumulation (MT merge pattern)

```cpp
// In ActionInitialization::Build(), RunAction::GenerateRun() returns a custom run:
G4Run* YourRunAction::GenerateRun() {
    fYourRun = new YourRun();
    return fYourRun;
}

// YourRun::Merge() is called by master to accumulate all worker runs:
void YourRun::Merge(const G4Run* run) {
    const YourRun* local = static_cast<const YourRun*>(run);
    fSum  += local->fSum;
    fSum2 += local->fSum2;
    G4Run::Merge(run);   // always call base last
}

// Access current run from EventAction:
YourRun* r = static_cast<YourRun*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());
```

### Messengers (UI command framework)

```cpp
class YourMessenger : public G4UImessenger {
public:
    YourMessenger(YourClass* obj) : fObj(obj) {
        fDir = new G4UIdirectory("/bbr/");
        fDir->SetGuidance("BBRsim commands");

        fCmd = new G4UIcmdWithADoubleAndUnit("/bbr/setZ", this);
        fCmd->SetParameterName("Z", false);   // false = not omittable
        fCmd->SetUnitCategory("Length");
        fCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
        fCmd->SetToBeBroadcasted(false);      // geometry: workers don't need it
    }
    ~YourMessenger() { delete fCmd; delete fDir; }

    void SetNewValue(G4UIcommand* cmd, G4String val) override {
        if (cmd == fCmd)
            fObj->SetZ(fCmd->GetNewDoubleValue(val));
    }
private:
    YourClass* fObj;
    G4UIdirectory*              fDir;
    G4UIcmdWithADoubleAndUnit*  fCmd;
};
```

Common command types: `G4UIcmdWithADoubleAndUnit`, `G4UIcmdWithADouble`,
`G4UIcmdWithAnInteger`, `G4UIcmdWithAString`, `G4UIcmdWithABool`,
`G4UIcmdWithoutParameter`.

`SetToBeBroadcasted(false)` — geometry-changing commands don't need to go to
workers because `Construct()` is called on the master. Physics/scoring
commands that affect worker-local objects should be `true` (default).

### Units and constants

```cpp
#include "G4SystemOfUnits.hh"   // mm, cm, m, eV, keV, MeV, GeV, ns, ps, deg, rad, K, ...
#include "G4PhysicalConstants.hh"  // c_light, h_Planck, k_Boltzmann, ...

// CLHEP equivalents (always valid):
#include "CLHEP/Units/SystemOfUnits.hh"
using CLHEP::mm; using CLHEP::MeV; // etc.

// Photon energy ↔ frequency: E = h*nu
G4double freq  = 500.e9;   // Hz
G4double energy = CLHEP::h_Planck * freq;   // in Geant4 internal units (MeV)
// Or: energy = 6.62607e-34 * freq / 1.602e-13 * MeV;
```

Geant4 internal length unit is **mm**; energy unit is **MeV**; time is **ns**.
Always multiply by unit when assigning: `G4double x = 3.0*cm;`.

### Key optical photon facts

- Particle name: `"opticalphoton"` — retrieved via
  `G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton")`.
- `G4ParticleGun` for optical photons: set `SetParticleEnergy(E)` where E
  is the photon energy (not kinetic energy in the usual sense).
- `SetParticleMomentumDirection` sets the k-direction; `SetParticlePolarization`
  sets the E-field polarization vector (must be perpendicular to k).
- `G4OpticsPhysics` (or `G4OpticalPhysics`) registers all optical processes:
  `G4OpAbsorption`, `G4OpRayleigh`, `G4OpMieHG`, `G4OpBoundaryProcess`,
  `G4OpWLS`, `G4Scintillation`, `G4Cerenkov`.
- To disable all optical processes except boundary (BBRsim ray-trace mode):
  use `G4OpticalPhysics` then call `optPhys->Configure(kAbsorption, false)`,
  etc., or disable them in the physics constructor.
- `RINDEX` must be defined for **both** materials at a boundary for Fresnel
  to work. If one material has no `RINDEX`, the photon is absorbed at the boundary.

### Touchable / volume navigation

```cpp
// In stepping action — find volume name:
G4String volName = step->GetPreStepPoint()
    ->GetTouchableHandle()->GetVolume()->GetName();

// Get material:
G4Material* mat = step->GetPreStepPoint()->GetMaterial();

// Get volume copy number:
G4int copyNo = step->GetPreStepPoint()
    ->GetTouchableHandle()->GetCopyNumber();

// Get rotation of the volume in world frame (used in BBRsim for crack orientation):
G4ThreeVector localDir = step->GetPreStepPoint()
    ->GetTouchableHandle()->GetHistory()->GetTopTransform()
    .TransformAxis(worldDir);
```

### Geometry reinitialisation (messenger-driven changes)

When a messenger command changes detector geometry at runtime:
```
/run/reinitializeGeometry   # full re-construction (needed when Construct() recomputes sizes)
/run/geometryModified        # lighter: flags geometry as dirty without re-running Construct()
```

In `SetTargetThickness()`, store the new value; the next `Initialize()` call
will re-run `Construct()` which reads the stored value.

### Debugging tips

- `G4cout` / `G4cerr` — use these instead of `std::cout`; they integrate with
  the Geant4 verbosity system.
- `/run/verbose 2`, `/event/verbose 2`, `/tracking/verbose 2` — progressively
  more stepping output in macro files.
- `/process/list` — lists all registered processes for all particles.
- Fixed-seed reproducibility: `G4Random::setTheSeed(seed)` before
  `BeamOn()`; or in macro: `/random/setSeeds seed1 seed2`.
