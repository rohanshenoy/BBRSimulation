# Next Steps Roadmap

> **Not an implementation plan** — a prioritised backlog of upcoming work. Convert individual items to full plans when ready to execute.

**Design constraint:** Do NOT patch `G4OpBoundaryProcess`. If a photon falls through the `BBSimOpBoundaryProcess` wrapper (i.e. the surface is not `vacuum_wg` and has no `REFLECTIVITY` on its MPT), let the stock Fresnel boundary process handle it. This is correct behaviour for dielectric boundaries (Cirlex→Si, vacuum→Cirlex, etc.).

---

## 1. Cirlex material  *(unblocked, data in hand)*

Add `GetCirlex()` to `include/BBRMaterials.hh`.

**Data source:** `build/yyc_cirlex_optical.csv` — digitised from YYC's `CDMSSnolabIRBackgroundMaterial.cc` (CDMSSnolabIRBackground, Zenbook era). 56-point tables on the YYC frequency grid (1 MHz – 100 THz):
- `n_Cirlex` — real refractive index, ~1.840 flat below 400 GHz, structured above
- `abslength_cm` — spans ~1e6 cm at MHz to ~0.01 cm at 10 THz

Material base: `G4_KAPTON` (YYC's approximation for Cirlex). Stock `G4OpBoundaryProcess` handles vacuum→Cirlex Fresnel automatically once RINDEX is set.

---

## 2. Si and Ge materials  *(need absorption length data)*

Add `GetSiliconCrystal()` and `GetGermaniumCrystal()` to `BBRMaterials.hh`.

**RINDEX:** Si=3.39, Ge≈4.0 (flat in trans-mm band, from YYC / Frey+ NASA Goddard).

**ABSLENGTH:** Derive from loss tangent (Chang §5.3.1.4):
- Si: tan δ = 1e-4  →  α = 2π n tan δ / λ  →  ABSLENGTH = 1/α
- Ge: tan δ = 6e-5

These are the detector materials. Correct ABSLENGTH is required for leakage-current post-processing (photon termination events in Si/Ge weighted by shallow-impurity MFP vs loss-tangent MFP).

---

## 3. Validation geometry  *(blocked — awaiting TAMU DR geometry)*

O2 goal: model the temperature-controlled OFHC-Cu BB source box + mesh-TES detector geometry at TAMU's DR.

Needs from Sunil/Nader:
- CAD or sketch of the BB source enclosure (dimensions, aperture, temperature range 4–20 K)
- Detector placement and active area
- Whether to use CADMesh (.STL import) or hand-coded G4 geometry for the first pass

Once geometry is defined: write `BBRSourceDetectorConstruction`, a validation mac, and ROOT TTree output for leakage-current post-processing.

---

## Already done (do not re-implement)

- Pass-through wrapper: `BBSimPhysics` + `BBSimOpBoundaryProcess`
- HFSS diffraction: `BBRHFSSData`, `HandleDiffractionBoundary`, `BBRCrackLibrary`
- Cu reflectance: `HandleReflectanceBoundary`, `GetOFHCCopper`, `GetOFCopperSerov`, `GetHPCopperSerov`
- Planck emitter: `BBRTestPGA`, `GetBBSpecCDF`, `ThermalSurface`, `BBEvt`
- Perfect reflector / absorber in `BBRMaterials.hh`
- YYC reference data extracted to `build/yyc_brass_reflectivity.csv`, `build/yyc_cirlex_optical.csv`
