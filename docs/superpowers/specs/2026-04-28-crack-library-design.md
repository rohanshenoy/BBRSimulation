# BBRCrackLibrary Design Spec
**Date:** 2026-04-28
**Status:** Approved

## Goal

Extract HFSS dataset routing out of `BBSimOpBoundaryProcess` into a dedicated
`BBRCrackLibrary` singleton, so the dispatch chain from photon step → correct
HFSS data is clean, documented, and generalizes to N cracks without any
per-crack registration.

## Context

Two crack volumes are now present in `BBRCrackDetectorConstruction`:

| Volume name                        | Position (world) | Gap b    | Depth   |
|------------------------------------|-----------------|----------|---------|
| `InfParallelPlate_crack1Rohan_500GHz` | (0, 0, 0) mm   | 50 µm   | 1.0 mm  |
| `InfParallelPlate_crack2_500GHz`      | (0, 0, 3) mm   | 100 µm  | 1.5 mm  |

Both volumes use material `vacuum_wg`. The boundary process fires on any
`vacuum_wg` entry; the volume name is the sole routing key.

No `G4OpticalSurface` / `G4LogicalBorderSurface` objects are used. The
material flag generalizes naturally to N cracks placed anywhere.

## Dispatch Chain

```
photon step
  └─ BBSimOpBoundaryProcess::PostStepDoIt
         └─ post-step material == "vacuum_wg"?
                └─ HandleDiffractionBoundary
                       └─ datasetId = volName.substr(0, volName.find(':'))
                              └─ BBRCrackLibrary::Instance().Lookup(datasetId)
                                     └─ lazy-load BBRHFSSData if not cached
                                            └─ transmittance + direction + position sampling
```

`BBRCrackLibrary` owns all I/O and caching. `BBSimOpBoundaryProcess` owns
all physics. `BBRHFSSData` is unchanged.

## BBRCrackLibrary Interface

```cpp
class BBRCrackLibrary {
 public:
  static BBRCrackLibrary& Instance();

  void         SetDataDir(const G4String& dir);   // default: "../HFSSSimData"
  BBRHFSSData& Lookup(const G4String& datasetId); // lazy-loads on first call

 private:
  BBRCrackLibrary() = default;
  G4String fDataDir = "../HFSSSimData";
  std::map<G4String, std::unique_ptr<BBRHFSSData>> fCache;
};
```

`Lookup` checks `fCache`; on miss constructs `BBRHFSSData(fDataDir, datasetId)`,
stores it, returns a reference. Adding crack N requires only placing a
`vacuum_wg` volume named after its dataset — no code changes.

`SetDataDir` allows the data path to be changed without recompiling
(e.g. from `BBRSim.cc` or a future messenger).

Thread safety: single-threaded only (`/run/numberOfThreads 1` enforced in
diffraction.mac). When MT is needed, `fCache` gets a mutex — change isolated
to this class.

## Changes to BBSimOpBoundaryProcess

`HandleDiffractionBoundary` replaces the inline cache block:

```cpp
// BEFORE
if (fHFSSCache.find(datasetId) == fHFSSCache.end())
    fHFSSCache[datasetId] = std::make_unique<BBRHFSSData>("../HFSSSimData", datasetId);
BBRHFSSData& hfss = *fHFSSCache[datasetId];

// AFTER
BBRHFSSData& hfss = BBRCrackLibrary::Instance().Lookup(datasetId);
```

`fHFSSCache` is removed from `BBSimOpBoundaryProcess.hh`. Everything else
in `HandleDiffractionBoundary` (angle decomposition, transmittance decision,
direction/position sampling) is untouched.

## Global vs. Local Coordinate Constraints

`BBRCrackLibrary` has no knowledge of coordinates. All frame handling stays
in `HandleDiffractionBoundary` using `xf = touch->GetHistory()->GetTopTransform()`
(maps world → local).

### Axis extraction
`normal_hat`, `theta_hat`, `phi_hat` come from `xf.InverseTransformAxis(local_axis)`,
which maps local → world (rotation only, translation-free). Correct for any
crack placement.

### Angle computation
`iwaveTheta` and `iwavePhi` are dot products of world-frame `khat` against the
world-frame expressions of the crack's local axes. Direction-only; translation
does not contaminate them. For crack2 at z=3 mm, angles must be computed
relative to crack2's local frame — this is guaranteed by `InverseTransformAxis`
extracting the correct rotation from `xf`.

### Position computation
`xf.TransformPoint(entryGl)` converts world entry point to local (subtracts
translation, applies inverse rotation). `DistanceToOut` operates in local
frame. `xf.InverseTransformPoint(localExit)` converts back to world (correctly
adds z=3 mm for crack2). `SampleExitPosition` uses world-frame `exitCenter`
and world-frame `theta_hat`/`phi_hat` — translation-independent. ✓

## Verification Requirements

1. **crack1 regression** — existing `diffraction.mac` (1000 events, z=0) must
   still give T_obs ≈ 0.528 after the refactor. Confirms the library dispatch
   produces identical results to the old inline cache.

2. **crack2 exit position** — run a photon at (−20, 0, 3 mm) along +x through
   crack2. Confirm `diffraction_output.csv` `pos_z_m` values cluster near
   3×10⁻³ m (world z ≈ 3 mm), not near 0. Confirms global/local position
   transforms are correct for the offset volume.

3. **crack2 oblique angle** — fire a photon at oblique incidence toward crack2.
   Confirm `iwaveTheta` / `iwavePhi` match the expected values from the crack's
   local geometry, not the world origin. Confirms angle computation is
   rotation-aware and translation-independent.

## Files Changed

| File | Action |
|------|--------|
| `include/BBRCrackLibrary.hh` | New |
| `src/BBRCrackLibrary.cc` | New |
| `src/BBSimOpBoundaryProcess.cc` | Remove inline cache, add `Lookup` call |
| `include/BBSimOpBoundaryProcess.hh` | Remove `fHFSSCache` member |

## Out of Scope

- Messenger for `SetDataDir` (future)
- Multi-threaded cache (future)
- Option C (volume-property tag routing) — deferred
- Angle/polarization sweep harness — separate feature
