# BBRCrackLibrary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract HFSS dataset routing from `BBSimOpBoundaryProcess` into a `BBRCrackLibrary` singleton so the dispatch chain is clean and generalizes to N cracks by volume-name matching.

**Architecture:** `BBRCrackLibrary` is a Meyer's singleton owning a lazy-loaded `std::map<datasetId, BBRHFSSData>` and a configurable base path. `BBSimOpBoundaryProcess::HandleDiffractionBoundary` calls `BBRCrackLibrary::Instance().Lookup(datasetId)` instead of managing its own cache. All coordinate math stays in `HandleDiffractionBoundary` unchanged.

**Tech Stack:** Geant4, C++17, CMake (`file(GLOB)` auto-picks up new `.cc` files — no CMakeLists change needed).

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `include/BBRCrackLibrary.hh` | Create | Singleton interface: `Instance`, `SetDataDir`, `Lookup` |
| `src/BBRCrackLibrary.cc` | Create | Lazy-load logic, cache, singleton definition |
| `include/BBSimOpBoundaryProcess.hh` | Modify | Remove `fHFSSCache`, remove `BBRHFSSData` include |
| `src/BBSimOpBoundaryProcess.cc` | Modify | Add `BBRCrackLibrary` include, replace 3-line cache block |
| `include/BBRDiffractionPGA.hh` | Modify | Add `gunZ_mm` constructor param (default 0) |
| `src/BBRDiffractionPGA.cc` | Modify | Use `gunZ_mm` in `SetParticlePosition` |
| `include/BBRDiffractionActionInit.hh` | Modify | Add `gunZ_mm` constructor param, `fGunZ_mm` member |
| `src/BBRDiffractionActionInit.cc` | Modify | Forward `fGunZ_mm` to `BBRDiffractionPGA` |
| `BBRSim.cc` | Modify | Detect `crack2` in mac name → pass `3.*mm` to ActionInit |
| `diffraction_crack2.mac` | Create | Smoke test firing at z=3 mm (crack2) |

---

## Task 1: Create `BBRCrackLibrary` header

**Files:**
- Create: `include/BBRCrackLibrary.hh`

- [ ] **Step 1: Write the header**

```cpp
#ifndef BBRCrackLibrary_hh
#define BBRCrackLibrary_hh

#include "BBRHFSSData.hh"
#include "G4String.hh"

#include <map>
#include <memory>

// Singleton that owns HFSS dataset loading and caching.
// Routing key: volume name stripped of any ":N" placement suffix = dataset ID.
// Lazy-loads BBRHFSSData on first Lookup; safe for single-threaded runs.
class BBRCrackLibrary
{
 public:
  static BBRCrackLibrary& Instance();

  void         SetDataDir(const G4String& dir);    // default: "../HFSSSimData"
  BBRHFSSData& Lookup(const G4String& datasetId);  // lazy-loads on first call

 private:
  BBRCrackLibrary() = default;
  G4String fDataDir = "../HFSSSimData";
  std::map<G4String, std::unique_ptr<BBRHFSSData>> fCache;
};

#endif
```

- [ ] **Step 2: Commit**

```bash
git add include/BBRCrackLibrary.hh
git commit -m "Add BBRCrackLibrary header — singleton HFSS dataset router"
```

---

## Task 2: Implement `BBRCrackLibrary`

**Files:**
- Create: `src/BBRCrackLibrary.cc`

- [ ] **Step 1: Write the implementation**

```cpp
#include "BBRCrackLibrary.hh"

BBRCrackLibrary& BBRCrackLibrary::Instance()
{
  static BBRCrackLibrary sInstance;
  return sInstance;
}

void BBRCrackLibrary::SetDataDir(const G4String& dir)
{
  fDataDir = dir;
}

BBRHFSSData& BBRCrackLibrary::Lookup(const G4String& datasetId)
{
  auto it = fCache.find(datasetId);
  if (it == fCache.end())
    it = fCache.emplace(datasetId,
                        std::make_unique<BBRHFSSData>(fDataDir, datasetId)).first;
  return *it->second;
}
```

- [ ] **Step 2: Commit**

```bash
git add src/BBRCrackLibrary.cc
git commit -m "Implement BBRCrackLibrary singleton — lazy HFSS dataset cache"
```

---

## Task 3: Wire `BBSimOpBoundaryProcess` to use `BBRCrackLibrary`

**Files:**
- Modify: `include/BBSimOpBoundaryProcess.hh`
- Modify: `src/BBSimOpBoundaryProcess.cc`

- [ ] **Step 1: Update the header — remove the cache member**

In `include/BBSimOpBoundaryProcess.hh`, make these changes:

Remove the `BBRHFSSData` include and the `<map>` / `<memory>` includes (no longer needed in the header):
```cpp
// REMOVE these three lines:
#include "BBRHFSSData.hh"
#include <map>
#include <memory>
```

Remove the `fHFSSCache` data member and its comment:
```cpp
// REMOVE:
  // Keyed by dataset ID (= volume name, stripped of any ":N" instance suffix).
  // Loaded lazily on first encounter; safe for single-threaded runs.
  std::map<G4String, std::unique_ptr<BBRHFSSData>> fHFSSCache;
```

The class comment at the top of the file should be updated to:
```cpp
// Wrapper around G4OpBoundaryProcess. PostStepDoIt intercepts photons that
// cross into a volume filled with material "vacuum_wg" and routes them through
// the HFSS diffraction model. The volume name is the HFSS dataset ID used to
// look up the correct CSV data via BBRCrackLibrary. Everything else is a
// pure pass-through.
```

- [ ] **Step 2: Update the .cc — replace inline cache with Lookup**

In `src/BBSimOpBoundaryProcess.cc`, add the include near the top (after existing includes):
```cpp
#include "BBRCrackLibrary.hh"
```

Find and replace the lazy-load block in `HandleDiffractionBoundary`:
```cpp
// BEFORE (3 lines to remove):
  if (fHFSSCache.find(datasetId) == fHFSSCache.end())
    fHFSSCache[datasetId] = std::make_unique<BBRHFSSData>("../HFSSSimData", datasetId);
  BBRHFSSData& hfss = *fHFSSCache[datasetId];

// AFTER (1 line):
  BBRHFSSData& hfss = BBRCrackLibrary::Instance().Lookup(datasetId);
```

- [ ] **Step 3: Commit**

```bash
git add include/BBSimOpBoundaryProcess.hh src/BBSimOpBoundaryProcess.cc
git commit -m "Wire BBSimOpBoundaryProcess to BBRCrackLibrary — remove inline cache"
```

---

## Task 4: Build and run crack1 regression

**Files:** (build only)

- [ ] **Step 1: Build**

```bash
cd build && make -j4 2>&1 | tail -5
```

Expected: `[100%] Linking CXX executable BBRSim` with no errors.

- [ ] **Step 2: Run crack1 smoke test**

```bash
cd build && ./BBRSim diffraction.mac 2>&1 | grep "T_obs"
```

Expected output (every 100 events):
```
[BBR] diffraction events=100 T_obs=0.5...
...
[BBR] diffraction events=1000 T_obs=0.5...
```

Final T_obs must be between 0.48 and 0.58 (theory ≈ 0.528).

- [ ] **Step 3: Commit if clean**

```bash
git add -p   # nothing to stage — build-only verification
```

No commit needed; build verification is sufficient.

---

## Task 5: Parameterize PGA and ActionInit for crack2

**Files:**
- Modify: `include/BBRDiffractionPGA.hh`
- Modify: `src/BBRDiffractionPGA.cc`
- Modify: `include/BBRDiffractionActionInit.hh`
- Modify: `src/BBRDiffractionActionInit.cc`

- [ ] **Step 1: Add `gunZ_mm` to `BBRDiffractionPGA`**

`include/BBRDiffractionPGA.hh` — replace default constructor with:
```cpp
  explicit BBRDiffractionPGA(G4double gunZ_mm = 0.);
```

`src/BBRDiffractionPGA.cc` — change constructor signature and position line:
```cpp
// BEFORE:
BBRDiffractionPGA::BBRDiffractionPGA()
{
  fGun = new G4ParticleGun(1);
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
  fGun->SetParticlePosition(G4ThreeVector(-20.*mm, 0., 0.));

// AFTER:
BBRDiffractionPGA::BBRDiffractionPGA(G4double gunZ_mm)
{
  fGun = new G4ParticleGun(1);
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
  fGun->SetParticlePosition(G4ThreeVector(-20.*mm, 0., gunZ_mm));
```

- [ ] **Step 2: Add `gunZ_mm` to `BBRDiffractionActionInit`**

`include/BBRDiffractionActionInit.hh` — replace:
```cpp
class BBRDiffractionActionInit : public G4VUserActionInitialization
{
 public:
  BBRDiffractionActionInit()  = default;
  ~BBRDiffractionActionInit() override = default;

  void BuildForMaster() const override;
  void Build()          const override;
};
```
With:
```cpp
#include "G4Types.hh"

class BBRDiffractionActionInit : public G4VUserActionInitialization
{
 public:
  explicit BBRDiffractionActionInit(G4double gunZ_mm = 0.);
  ~BBRDiffractionActionInit() override = default;

  void BuildForMaster() const override;
  void Build()          const override;

 private:
  G4double fGunZ_mm = 0.;
};
```

`src/BBRDiffractionActionInit.cc` — add constructor and update Build:
```cpp
#include "BBRDiffractionActionInit.hh"

#include "BBRDiffractionPGA.hh"
#include "RunAction.hh"

BBRDiffractionActionInit::BBRDiffractionActionInit(G4double gunZ_mm)
  : fGunZ_mm(gunZ_mm) {}

void BBRDiffractionActionInit::BuildForMaster() const
{
  SetUserAction(new RunAction(nullptr));
}

void BBRDiffractionActionInit::Build() const
{
  SetUserAction(new BBRDiffractionPGA(fGunZ_mm));
  SetUserAction(new RunAction(nullptr));
}
```

- [ ] **Step 3: Commit**

```bash
git add include/BBRDiffractionPGA.hh src/BBRDiffractionPGA.cc \
        include/BBRDiffractionActionInit.hh src/BBRDiffractionActionInit.cc
git commit -m "Parameterize BBRDiffractionPGA/ActionInit with gunZ_mm for crack2 test"
```

---

## Task 6: Wire `BBRSim.cc` and add `diffraction_crack2.mac`

**Files:**
- Modify: `BBRSim.cc`
- Create: `diffraction_crack2.mac`

- [ ] **Step 1: Update `BBRSim.cc` to detect crack2 mac**

Add `#include "G4SystemOfUnits.hh"` to the includes block in `BBRSim.cc`.

Replace the ActionInit selection block:
```cpp
// BEFORE:
  if (argc > 1 && G4String(argv[1]).find("diffraction") != G4String::npos)
    runManager->SetUserInitialization(new BBRDiffractionActionInit());
  else
    runManager->SetUserInitialization(new ActionInitialization());

// AFTER:
  if (argc > 1 && G4String(argv[1]).find("diffraction") != G4String::npos) {
    G4double gunZ = (G4String(argv[1]).find("crack2") != G4String::npos)
                    ? 3.*mm : 0.;
    runManager->SetUserInitialization(new BBRDiffractionActionInit(gunZ));
  } else {
    runManager->SetUserInitialization(new ActionInitialization());
  }
```

- [ ] **Step 2: Create `diffraction_crack2.mac`**

```
# Diffraction smoke test — crack2 (b=100 µm, depth=1.5 mm).
# Gun at (-20, 0, 3) mm along +x — hits crack2 centered at (x=0, y=0, z=3 mm).
# Verification: diffraction_output.csv pos_z_m should cluster near 3e-3 m.

/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/run/numberOfThreads 1
/run/initialize

/random/setSeeds 1234 5678

/run/beamOn 1000
```

- [ ] **Step 3: Commit**

```bash
git add BBRSim.cc diffraction_crack2.mac
git commit -m "Wire BBRSim.cc for crack2 test mac; add diffraction_crack2.mac"
```

---

## Task 7: Build and verify crack2 global/local coordinates

**Files:** (build + run + verify)

- [ ] **Step 1: Build**

```bash
cd build && make -j4 2>&1 | tail -5
```

Expected: clean build, no errors.

- [ ] **Step 2: Run crack2 smoke test**

```bash
cd build && ./BBRSim diffraction_crack2.mac 2>&1 | grep "T_obs" | tail -3
```

Expected: T_obs some value printed (crack2 transmittance — may differ from 0.528). No crash.

- [ ] **Step 3: Verify exit positions cluster near z=3 mm**

```bash
conda run -n bbrsim python - <<'EOF'
import pandas as pd, numpy as np
df = pd.read_csv('build/diffraction_output.csv')
# Only rows from the crack2 run (rerun from scratch)
mean_z = df['pos_z_m'].mean()
std_z  = df['pos_z_m'].std()
print(f"pos_z_m  mean={mean_z:.4e}  std={std_z:.4e}")
assert abs(mean_z - 3e-3) < 5e-4, f"FAIL: expected ~3e-3, got {mean_z:.4e}"
print("PASS: exit positions are near crack2 center (z=3 mm)")
EOF
```

Expected output:
```
pos_z_m  mean=3.????e-03  std=?.????e-05
PASS: exit positions are near crack2 center (z=3 mm)
```

- [ ] **Step 4: Re-run crack1 regression to confirm nothing regressed**

```bash
cd build && ./BBRSim diffraction.mac 2>&1 | grep "T_obs" | tail -1
```

Expected: T_obs between 0.48 and 0.58.

- [ ] **Step 5: Final commit**

```bash
# All code already committed in prior tasks — tag the verification as clean
git log --oneline -6
```

No new commit needed unless any fixups were required during verification.

---

## Verification Checklist (from spec)

- [ ] crack1 regression: T_obs ≈ 0.528 after refactor ← Task 4
- [ ] crack2 exit positions: `pos_z_m` ≈ 3×10⁻³ m ← Task 7 Step 3
- [ ] crack2 oblique angle: manual inspection of `[BBR]` log lines at oblique incidence (future — not automated in this plan)

## Out of Scope

- Messenger for `SetDataDir`
- Multi-threaded cache
- Option C (volume-property tag routing)
- Oblique-angle automated test (noted in spec; manual verification deferred)
