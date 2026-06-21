#ifndef BBRTestDetectorConstruction_hh
#define BBRTestDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

// Test geometry: 50 cm world, 4 mm Cu slab at x=2 mm with two vacuum_wg crack
// daughters. The Cu wall material is built from BBRConfigManager's CuRRR /
// CuStageT_K at Construct() time via BBRMaterials::GetCopper(RRR, T_K).
//
// Configuration commands (owned by BBRConfigMessenger; before /run/initialize):
//   /bbr/det/setCuMaterial <OFHC_Cu|OF_Cu|HP_Cu>   named alias -> RRR
//   /bbr/det/setCuRRR <N>                           direct integer RRR (>= 1)
//   /bbr/det/setCuStageT <T> K                      temperature stage [K]
class BBRTestDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  BBRTestDetectorConstruction() = default;
  ~BBRTestDetectorConstruction() override = default;
  G4VPhysicalVolume* Construct() override;
};

#endif
