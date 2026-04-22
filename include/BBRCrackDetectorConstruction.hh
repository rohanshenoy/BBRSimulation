#ifndef BBRCrackDetectorConstruction_hh
#define BBRCrackDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

// Minimal geometry for the diffraction smoke test:
//   - Vacuum world (200 mm cube)
//   - TEM_waveguide_crack slab at origin matching HFSS model dimensions
//     (a=10 mm propagation depth, b=0.05 mm gap, long dimension truncated to 100 mm)
// Normal incidence photons travel along +x; BBSimOpBoundaryProcess intercepts
// the waveguide boundary and routes to HandleDiffractionBoundary.
class BBRCrackDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  BBRCrackDetectorConstruction()  = default;
  ~BBRCrackDetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;
};

#endif
