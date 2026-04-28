#ifndef BBRCrackDetectorConstruction_hh
#define BBRCrackDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

// Minimal geometry for the diffraction smoke test:
//   - Vacuum world (200 mm cube)
//   - crack1 (InfParallelPlate_crack1Rohan_500GHz) centered at (x=0, y=0, z=0):
//       HFSS ZSize=1 mm depth, XSize=0.05 mm gap, YSize=10 mm long
//   - crack2 (InfParallelPlate_crack2_500GHz) centered at (x=0, y=0, z=3 mm):
//       HFSS ZSize=1.5 mm depth, XSize=0.1 mm gap, YSize=10 mm long
// Both slabs sit at the same propagation depth (x=0), side by side in z.
// A photon at z=0 hits crack1; a photon offset to z=3 mm hits crack2.
// BBSimOpBoundaryProcess loads the per-crack HFSS dataset by volume name.
class BBRCrackDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  BBRCrackDetectorConstruction()  = default;
  ~BBRCrackDetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;
};

#endif
