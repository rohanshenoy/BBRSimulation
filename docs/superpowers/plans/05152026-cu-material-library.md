# Cu Material Library Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two experimentally-grounded copper materials (`OF_Cu` and `HP_Cu`) to `BBRMaterials.hh`, derived from Serov et al. (IEEE TMT 2016) mm-wave reflectivity data at 4 K, alongside a Python analytical pre-check and Geant4 integration.

**Architecture:** Both new materials use the same Hagen-Rubens formula already used by `OFHC_Cu` — `R = 1 − 2√(2ε₀ω/σ_eff)` — but with effective conductivities derived from Serov's 4 K experimental data rather than the RRR=100 theoretical value. `OF_Cu` uses σ_eff = 1/ρ₀ where ρ₀=0.56×10⁻⁸ Ω·m is Serov's measured residual resistance for 99.97% OF copper; `HP_Cu` back-calculates σ_eff from Serov's reported D=0.55×10⁻³ at 230 GHz for 99.999% hydrogen-annealed copper. `BBRTestDetectorConstruction` gains an optional `G4Material*` constructor argument; `BBRSim.cc` selects the material based on the mac filename. `check_cu_absorptance.py` is generalised to accept any material name; a new `check_cu_serov.py` validates the σ_eff values analytically against Serov's published data points before any C++ is touched.

**Physics reference:** Serov, Parshin, Bubnov, IEEE TMT 64(11), 3828 (2016). Key measurements at T=4 K:
- OF copper (99.97%): D ≈ 0.58×10⁻³ at 150 GHz → σ_eff = 1/ρ₀ = 1.786×10⁸ S/m
- HP copper annealed (99.999%): D ≈ 0.55×10⁻³ at 230 GHz → σ_eff = 3.38×10⁸ S/m (back-calc)
- Existing OFHC_Cu: σ = 5.96×10⁹ S/m (RRR=100 theoretical, optimistic lower bound on loss)

**Tech Stack:** Geant4 11.4, C++17, CLHEP, Python 3 (conda env `bbrsim`), numpy, scipy.

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `scripts/check_cu_serov.py` | Create | Analytical: verify σ_eff values reproduce Serov's D(f) at T=4K |
| `include/BBRMaterials.hh` | Modify | Add `GetOFCopperSerov()` and `GetHPCopperSerov()` |
| `include/BBRTestDetectorConstruction.hh` | Modify | Add `G4Material* fCuMat` member + constructor arg |
| `src/BBRTestDetectorConstruction.cc` | Modify | Use `fCuMat` instead of hard-coded `GetOFHCCopper()` |
| `BBRSim.cc` | Modify | Detect mac filename; pass correct Cu material to constructor |
| `test_of_cu.mac` | Create | 1 000-event smoke test with OF_Cu slab |
| `test_hp_cu.mac` | Create | 1 000-event smoke test with HP_Cu slab |
| `CMakeLists.txt` | Modify | Copy new `.mac` files to build dir |
| `scripts/check_cu_absorptance.py` | Modify | Generalise regex to accept any `mat=<name>` |

---

## Task 0: Analytical pre-check — write and pass check_cu_serov.py

Write a pure-Python script that verifies our σ_eff values reproduce Serov's 4 K measurements before any C++ is written. Run it and confirm PASS.

**Files:**
- Create: `scripts/check_cu_serov.py`

- [ ] **Step 1: Write `scripts/check_cu_serov.py`**

```python
#!/usr/bin/env python3
"""
check_cu_serov.py
Verify that the Hagen-Rubens σ_eff values chosen for OF_Cu and HP_Cu
reproduce Serov et al. (IEEE TMT 2016) reflection-loss measurements at T=4 K.

Reference points (Serov Figs 6 and 8, read at T=4 K):
  OF copper  (99.97%),          150 GHz: D = 0.58e-3
  HP copper  (99.999%, annealed), 230 GHz: D = 0.55e-3

Run: conda run -n bbrsim python scripts/check_cu_serov.py
"""
import numpy as np
import sys

eps0 = 8.8541878128e-12   # F/m

def hagen_rubens_D(freq_Hz, sigma_SI):
    omega = 2.0 * np.pi * freq_Hz
    return 2.0 * np.sqrt(2.0 * eps0 * omega / sigma_SI)

# ------------------------------------------------------------------
# 1.  OF copper: σ_eff = 1/ρ₀  (Serov fit: ρ₀=0.56e-8 Ω·m)
# ------------------------------------------------------------------
rho0_OF   = 0.56e-8          # Ω·m  (Serov Eq. 24 fit to OF copper)
sigma_OF  = 1.0 / rho0_OF   # S/m
D_OF_calc = hagen_rubens_D(150e9, sigma_OF)
D_OF_ref  = 0.58e-3          # Serov Fig. 6, T=4 K, 150 GHz

# ------------------------------------------------------------------
# 2.  HP copper annealed: σ_eff back-calculated from D=0.55e-3 at 230 GHz
#     σ = 8 ε₀ ω / D²  (inverted Hagen-Rubens)
# ------------------------------------------------------------------
D_HP_ref   = 0.55e-3         # Serov Fig. 8, T=4 K, 230 GHz
omega_230  = 2.0 * np.pi * 230e9
sigma_HP   = 8.0 * eps0 * omega_230 / D_HP_ref**2
D_HP_calc  = hagen_rubens_D(230e9, sigma_HP)

# ------------------------------------------------------------------
# 3.  Existing OFHC_Cu (theoretical, for reference)
# ------------------------------------------------------------------
sigma_OFHC = 5.96e9
D_OFHC_500 = hagen_rubens_D(500e9, sigma_OFHC)

print("=" * 60)
print(f"{'Material':<20} {'freq':>8} {'D_calc':>10} {'D_ref':>10} {'ratio':>8}")
print("-" * 60)

ratio_OF = D_OF_calc / D_OF_ref
ratio_HP = D_HP_calc / D_HP_ref
tol = 0.10   # ±10 % tolerance

rows = [
    ("OF_Cu (Serov)",  "150 GHz", D_OF_calc, D_OF_ref, ratio_OF),
    ("HP_Cu (Serov)",  "230 GHz", D_HP_calc, D_HP_ref, ratio_HP),
    ("OFHC_Cu (ref)", "500 GHz", D_OFHC_500, None, None),
]
for name, freq, calc, ref, ratio in rows:
    ref_str   = f"{ref:.3e}" if ref  is not None else "   N/A   "
    ratio_str = f"{ratio:.3f}" if ratio is not None else "  N/A  "
    print(f"{name:<20} {freq:>8} {calc:>10.3e} {ref_str:>10} {ratio_str:>8}")

print("=" * 60)
print(f"\nσ_eff  OF_Cu  = {sigma_OF:.4e} S/m   (= 1/ρ₀, ρ₀={rho0_OF:.2e} Ω·m)")
print(f"σ_eff  HP_Cu  = {sigma_HP:.4e} S/m   (back-calc, D=0.55e-3 at 230 GHz)")
print(f"σ      OFHC_Cu= {sigma_OFHC:.4e} S/m   (RRR=100 × σ_RT, theoretical)")

passed = (abs(ratio_OF - 1.0) <= tol) and (abs(ratio_HP - 1.0) <= tol)
print(f"\nRESULT: {'PASS' if passed else 'FAIL'}  (tolerance ±{int(tol*100)}%)")
sys.exit(0 if passed else 1)
```

- [ ] **Step 2: Run the script**

```bash
conda run -n bbrsim python scripts/check_cu_serov.py
```

Expected output (both ratios within ±10%):
```
====================================================
Material              freq     D_calc      D_ref    ratio
----------------------------------------------------
OF_Cu (Serov)      150 GHz  6.110e-04  5.800e-04    1.053
HP_Cu (Serov)      230 GHz  5.500e-04  5.500e-04    1.000
OFHC_Cu (ref)      500 GHz  1.932e-04       N/A      N/A
====================================================
σ_eff  OF_Cu  = 1.7857e+08 S/m
σ_eff  HP_Cu  = 3.3832e+08 S/m
σ      OFHC_Cu= 5.9600e+09 S/m
RESULT: PASS  (tolerance ±10%)
```

- [ ] **Step 3: Commit**

```bash
git add scripts/check_cu_serov.py
git commit -m "feat: analytical pre-check for OF_Cu and HP_Cu sigma_eff vs Serov (2016)"
```

---

## Task 1: Add GetOFCopperSerov() and GetHPCopperSerov() to BBRMaterials.hh

**Files:**
- Modify: `include/BBRMaterials.hh`

- [ ] **Step 1: Add the two new inline factory functions**

Open `include/BBRMaterials.hh`. After the closing `}` of `GetPerfectReflector()` and before the closing `} // namespace BBRMaterials`, insert:

```cpp
// OF copper at 4 K (Serov et al., IEEE TMT 2016, Fig. 6).
// 99.97% Cu, residual-resistance dominated: σ_eff = 1/ρ₀, ρ₀=0.56e-8 Ω·m.
// Hagen-Rubens D(f) ≈ 0.61e-3 at 150 GHz (Serov: 0.58e-3, within 6%).
inline G4Material* GetOFCopperSerov()
{
  G4Material* mat = G4Material::GetMaterial("OF_Cu", false);
  if (mat) return mat;
  mat = new G4Material("OF_Cu", 29., 63.546*g/mole, 8.96*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;        // 50 GHz
  const G4double Emax  = 8.27e-2*eV;        // 20 THz
  const G4double eps0  = 8.8541878128e-12;  // F/m
  const G4double sigma = 1.786e8;           // S/m  (= 1/ρ₀, ρ₀=0.56e-8 Ω·m)
  const G4double h_eVs = 4.13566769692e-15; // eV·s

  G4double energies[N], refls[N];
  G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma);
    refls[i]       = std::max(0., std::min(1., R));
  }
  G4double e2[]     = {Emin, Emax};
  G4double rindex[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       rindex, 2);
  mpt->AddProperty("REFLECTIVITY", energies, refls,  N);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// HP copper (hydrogen-annealed) at 4 K (Serov et al., IEEE TMT 2016, Fig. 8).
// 99.999% Cu, annealed; σ_eff back-calculated from D=0.55e-3 at 230 GHz.
// Hagen-Rubens D(f) ≈ 0.55e-3 at 230 GHz by construction.
inline G4Material* GetHPCopperSerov()
{
  G4Material* mat = G4Material::GetMaterial("HP_Cu", false);
  if (mat) return mat;
  mat = new G4Material("HP_Cu", 29., 63.546*g/mole, 8.96*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;        // 50 GHz
  const G4double Emax  = 8.27e-2*eV;        // 20 THz
  const G4double eps0  = 8.8541878128e-12;  // F/m
  const G4double sigma = 3.383e8;           // S/m  (back-calc, D=0.55e-3 @ 230 GHz)
  const G4double h_eVs = 4.13566769692e-15; // eV·s

  G4double energies[N], refls[N];
  G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma);
    refls[i]       = std::max(0., std::min(1., R));
  }
  G4double e2[]     = {Emin, Emax};
  G4double rindex[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       rindex, 2);
  mpt->AddProperty("REFLECTIVITY", energies, refls,  N);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}
```

- [ ] **Step 2: Verify the file compiles (header-only, so build the project)**

```bash
cd build && cmake .. -DWITH_GEANT4_UIVIS=OFF && make -j$(sysctl -n hw.logicalcpu) 2>&1 | tail -20
```

Expected: no errors, `BBRSim` rebuilt.

- [ ] **Step 3: Commit**

```bash
cd ..
git add include/BBRMaterials.hh
git commit -m "feat: add OF_Cu and HP_Cu materials (Serov 2016 experimental conductivities at 4K)"
```

---

## Task 2: Add Cu material selector to BBRTestDetectorConstruction

Allow passing the Cu wall material at construction time. Default remains `GetOFHCCopper()` so all existing macs and tests are unaffected.

**Files:**
- Modify: `include/BBRTestDetectorConstruction.hh`
- Modify: `src/BBRTestDetectorConstruction.cc`

- [ ] **Step 1: Update the header**

Replace the contents of `include/BBRTestDetectorConstruction.hh` with:

```cpp
#ifndef BBRTestDetectorConstruction_hh
#define BBRTestDetectorConstruction_hh
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"

class BBRTestDetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit BBRTestDetectorConstruction(G4Material* cuMat = nullptr);
  ~BBRTestDetectorConstruction() override = default;
  G4VPhysicalVolume* Construct() override;
private:
  G4Material* fCuMat;  // nullptr → GetOFHCCopper() at construct time
};
#endif
```

- [ ] **Step 2: Update the source**

In `src/BBRTestDetectorConstruction.cc`, replace the first two lines of the file (the `#include`s before the function) plus the function opening with the version below. The key change: add the constructor definition and replace the hard-coded `BBRMaterials::GetOFHCCopper()` call with `fCuMat`.

The full file after editing:

```cpp
#include "BBRTestDetectorConstruction.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

BBRTestDetectorConstruction::BBRTestDetectorConstruction(G4Material* cuMat)
  : G4VUserDetectorConstruction()
  , fCuMat(cuMat ? cuMat : BBRMaterials::GetOFHCCopper())
{}

G4VPhysicalVolume* BBRTestDetectorConstruction::Construct()
{
  // World: 50 cm cube of G4_Galactic with RINDEX=1
  G4Material* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto* worldMPT = new G4MaterialPropertiesTable();
  G4double e[]  = {1e-6*eV, 1.0*eV};
  G4double ri[] = {1., 1.};
  worldMPT->AddProperty("RINDEX", e, ri, 2);
  vac->SetMaterialPropertiesTable(worldMPT);

  auto* worldBox = new G4Box("World", 250.*mm, 250.*mm, 250.*mm);
  auto* worldLV  = new G4LogicalVolume(worldBox, vac, "World");

  // Cu wall: center at (2mm,0,0), front face at x=0, back face at x=4mm
  auto* cuBox = new G4Box("CuSlab", 2.*mm, 25.*mm, 25.*mm);
  auto* cuLV  = new G4LogicalVolume(cuBox, fCuMat, "CuSlab");
  new G4PVPlacement(nullptr, G4ThreeVector(2.*mm, 0., 0.),
                    cuLV, "CuSlab", worldLV, false, 0, true);

  // crack1: full-span slab daughter so HFSS-transmitted photons exit into world at x=4mm.
  {
    static const char* kId = "InfParallelPlate_crack1Rohan_500GHz";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.026*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
                      logical, kId, cuLV, false, 0, true);
  }

  // crack2: full-span slab daughter.
  {
    static const char* kId = "InfParallelPlate_crack2_500GHz";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.051*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 3.*mm),
                      logical, kId, cuLV, false, 0, true);
  }

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLV, "World", nullptr, false, 0, true);
}
```

- [ ] **Step 3: Build**

```bash
cd build && make -j$(sysctl -n hw.logicalcpu) 2>&1 | tail -20
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
cd ..
git add include/BBRTestDetectorConstruction.hh src/BBRTestDetectorConstruction.cc
git commit -m "feat: BBRTestDetectorConstruction accepts optional Cu material arg (default OFHC_Cu)"
```

---

## Task 3: Route mac names to Cu material in BBRSim.cc

When the mac file contains `of_cu`, use `GetOFCopperSerov()`; when it contains `hp_cu`, use `GetHPCopperSerov()`; otherwise use the existing default.

**Files:**
- Modify: `BBRSim.cc`

- [ ] **Step 1: Add material-selection logic**

Replace the existing `BBRSim.cc` with:

```cpp
/// \file BBRSim.cc
/// \brief Main program for BBRSimulation. Adapted from Geant4 OpNovice2.

#include "BBRTestDetectorConstruction.hh"
#include "BBRTestActionInit.hh"
#include "BBRMaterials.hh"
#include "BBSimPhysics.hh"
#include "SteppingVerbose.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include <string>

int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  auto* steppingVerbose = new SteppingVerbose;
  auto* runManager      = G4RunManagerFactory::CreateRunManager();

  // Select Cu wall material from mac filename.
  G4Material* cuMat = nullptr;
  if (argc > 1) {
    std::string mac(argv[1]);
    if (mac.find("of_cu") != std::string::npos)
      cuMat = BBRMaterials::GetOFCopperSerov();
    else if (mac.find("hp_cu") != std::string::npos)
      cuMat = BBRMaterials::GetHPCopperSerov();
    // else nullptr → BBRTestDetectorConstruction defaults to OFHC_Cu
  }

  runManager->SetUserInitialization(new BBRTestDetectorConstruction(cuMat));

  auto* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->RegisterPhysics(new BBSimPhysics());
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new BBRTestActionInit());

  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (ui) {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    UImanager->ApplyCommand(G4String("/control/execute ") + G4String(argv[1]));
  }

  delete visManager;
  delete runManager;
  delete steppingVerbose;
  return 0;
}
```

- [ ] **Step 2: Build**

```bash
cd build && make -j$(sysctl -n hw.logicalcpu) 2>&1 | tail -20
```

Expected: no errors.

- [ ] **Step 3: Commit**

```bash
cd ..
git add BBRSim.cc
git commit -m "feat: BBRSim routes of_cu/hp_cu mac names to Serov Cu materials"
```

---

## Task 4: Add test macs and update CMakeLists

**Files:**
- Create: `test_of_cu.mac`
- Create: `test_hp_cu.mac`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Write `test_of_cu.mac`**

```mac
# test_of_cu.mac
# Smoke test for OF_Cu (Serov 2016 OF copper at 4 K, sigma_eff = 1/rho_0).
# Expect [BBR] reflectance A_obs ~ 1e-3 (higher loss than OFHC_Cu).
/run/numberOfThreads 1
/random/setSeeds 12345 67890
/bbr/gun/mode false
/bbr/thermal/setT 4.0
/run/initialize
/run/beamOn 1000
```

- [ ] **Step 2: Write `test_hp_cu.mac`**

```mac
# test_hp_cu.mac
# Smoke test for HP_Cu (Serov 2016 HP annealed copper at 4 K, sigma_eff back-calc).
# Expect [BBR] reflectance A_obs ~ 8e-4 (better than OF_Cu, worse than OFHC_Cu).
/run/numberOfThreads 1
/random/setSeeds 12345 67890
/bbr/gun/mode false
/bbr/thermal/setT 4.0
/run/initialize
/run/beamOn 1000
```

- [ ] **Step 3: Add the new macs to CMakeLists.txt**

Find the block in `CMakeLists.txt` that copies mac files to the build directory. It will look like a `configure_file` or `file(COPY ...)` call listing individual macs, or a GLOB. Add `test_of_cu.mac` and `test_hp_cu.mac` to the same list.

The copy block currently contains entries like `test_50M.mac` and `test_10K.mac`. Add:
```cmake
configure_file(test_of_cu.mac  ${CMAKE_CURRENT_BINARY_DIR}/test_of_cu.mac  COPYONLY)
configure_file(test_hp_cu.mac  ${CMAKE_CURRENT_BINARY_DIR}/test_hp_cu.mac  COPYONLY)
```

(Match the exact pattern already used in the file — if it uses `file(COPY ...)` with a list, add the two mac names to that list instead.)

- [ ] **Step 4: Build, then run both smoke tests**

```bash
cd build && cmake .. -DWITH_GEANT4_UIVIS=OFF && make -j$(sysctl -n hw.logicalcpu)
./BBRSim test_of_cu.mac 2>&1 | grep "\[BBR\]"
./BBRSim test_hp_cu.mac 2>&1 | grep "\[BBR\]"
```

Expected: each run prints `[BBR] reflectance mat=OF_Cu ...` / `mat=HP_Cu ...` lines. `A_obs` for `OF_Cu` should be higher than for `HP_Cu`, and both should be higher than the `A_obs` seen with `OFHC_Cu` in the default test.

- [ ] **Step 5: Commit**

```bash
cd ..
git add test_of_cu.mac test_hp_cu.mac CMakeLists.txt
git commit -m "feat: add test_of_cu.mac and test_hp_cu.mac smoke tests; update CMakeLists"
```

---

## Task 5: Generalise check_cu_absorptance.py and add per-material check

The current script hard-codes `mat=OFHC_Cu` in its regex. Generalise it so it works with any material name, and add `--sigma` so callers can pass the right σ_eff for the Planck-weighted theory.

**Files:**
- Modify: `scripts/check_cu_absorptance.py`

- [ ] **Step 1: Update the script**

Replace `scripts/check_cu_absorptance.py` with:

```python
"""
check_cu_absorptance.py
Parse BBRSim stdout for [BBR] reflectance lines; compare A_obs to
Hagen-Rubens Planck-weighted theory. PASS if 0.3 < A_obs/A_theory < 3.0.

Usage:
    ./BBRSim test.mac >bbrsim_out.txt 2>&1
    conda run -n bbrsim python scripts/check_cu_absorptance.py bbrsim_out.txt [--sigma S]

Options:
  --sigma S   Effective conductivity in S/m used to compute A_theory.
              Defaults to 5.96e9 (OFHC_Cu RRR=100).
              Use 1.786e8 for OF_Cu, 3.383e8 for HP_Cu.
"""
import sys
import re
import argparse
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="?", default="bbrsim_stdout.txt")
parser.add_argument("--sigma", type=float, default=5.96e9,
                    help="σ_eff in S/m for Hagen-Rubens theory (default: 5.96e9 = OFHC_Cu)")
args = parser.parse_args()

pattern = re.compile(
    r"\[BBR\] reflectance mat=(\S+) N=(\d+) A_obs=([\d.e+-]+) R_theory=([\d.e+-]+)"
)
best = None
with open(args.file) as f:
    for line in f:
        m = pattern.search(line)
        if m:
            mat, n, a_obs = m.group(1), int(m.group(2)), float(m.group(3))
            if best is None or n > best[1]:
                best = (mat, n, a_obs)

if best is None:
    print("ERROR: no [BBR] reflectance lines found in", args.file)
    sys.exit(1)

mat_name, N, A_obs = best
print(f"[BBR] log  mat={mat_name}  N={N}  A_obs={A_obs:.6e}")

eps0  = 8.854e-12
sigma = args.sigma
hbar  = 1.0546e-34
k_B   = 1.3806e-23
T     = 4.0

def A_HR(omega):
    return 2.0 * np.sqrt(2.0 * eps0 * omega / sigma)

def planck_weight(omega):
    x = hbar * omega / (k_B * T)
    return np.where(x < 500, omega**2 / (np.expm1(x)), 0.0)

omega_lo = 2 * np.pi * 1e9
omega_hi = 2 * np.pi * 20e12

num, _ = integrate.quad(lambda w: A_HR(w) * planck_weight(w), omega_lo, omega_hi)
den, _ = integrate.quad(lambda w:            planck_weight(w), omega_lo, omega_hi)
A_theory = num / den

ratio = A_obs / A_theory
lo, hi = 0.3, 3.0
passed = lo < ratio < hi

print(f"σ_eff used         = {sigma:.4e} S/m")
print(f"A_theory (Planck)  = {A_theory:.6e}")
print(f"A_obs              = {A_obs:.6e}")
print(f"ratio A_obs/theory = {ratio:.3f}  (expected [{lo}, {hi}])")
print(f"RESULT: {'PASS' if passed else 'FAIL'}")

sys.exit(0 if passed else 1)
```

- [ ] **Step 2: Run against the OF_Cu and HP_Cu smoke-test outputs**

```bash
cd build
./BBRSim test_of_cu.mac >of_cu_out.txt 2>&1
./BBRSim test_hp_cu.mac >hp_cu_out.txt 2>&1

conda run -n bbrsim python ../scripts/check_cu_absorptance.py of_cu_out.txt --sigma 1.786e8
conda run -n bbrsim python ../scripts/check_cu_absorptance.py hp_cu_out.txt --sigma 3.383e8
```

Expected: both print `RESULT: PASS`.

Also confirm the existing default test still passes:

```bash
./BBRSim test.mac >test_out.txt 2>&1
conda run -n bbrsim python ../scripts/check_cu_absorptance.py test_out.txt
```

Expected: `RESULT: PASS` (default σ=5.96e9 matches OFHC_Cu).

- [ ] **Step 3: Commit**

```bash
cd ..
git add scripts/check_cu_absorptance.py
git commit -m "feat: generalise check_cu_absorptance for any Cu material via --sigma flag"
```

---

## Self-Review

**Spec coverage:** All three materials defined. Analytical pre-check written. C++ materials added to BBRMaterials.hh with same pattern as existing. Constructor arg plumbed through. BBRSim.cc routes by mac name. Test macs created. CMakeLists updated. check_cu_absorptance.py generalised and validated against all three materials.

**Placeholder scan:** No TBDs. All code blocks are complete.

**Type consistency:** `G4Material*` parameter in constructor header matches usage in `.cc`. `GetOFCopperSerov()` / `GetHPCopperSerov()` names match between BBRMaterials.hh and BBRSim.cc. `--sigma` arg name matches across check_cu_absorptance.py steps.

**σ_eff values used consistently:**
- `OF_Cu`: `sigma = 1.786e8` in BBRMaterials.hh Task 1; `--sigma 1.786e8` in Task 5; `sigma_OF = 1/rho0_OF` in Task 0 (same value)
- `HP_Cu`: `sigma = 3.383e8` in BBRMaterials.hh Task 1; `--sigma 3.383e8` in Task 5; `sigma_HP` back-calc in Task 0 (same value by construction)
