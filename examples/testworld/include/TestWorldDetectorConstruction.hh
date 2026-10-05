#ifndef TestWorldDetectorConstruction_hh
#define TestWorldDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class TestWorldMessenger;

// Test geometry: 50 cm world, 4 mm Cu slab at x=2 mm with two vacuum_wg crack
// daughters, plus, only with /bbr/testworld/roundGap true, a Cu plate
// RoundGapPlate holding the straight round gap RoundGap_r50um (a vacuum_wg
// displaced tube). The Cu material is built from BBRConfigManager's CuRRR /
// CuStageT_K at Construct() time via BBRMaterials::GetCopper(RRR, T_K).
//
// Configuration commands (before /run/initialize):
//   from BBRConfigMessenger (library):
//     /bbr/det/setCuMaterial <OFHC_Cu|OF_Cu|HP_Cu>   named alias -> RRR
//     /bbr/det/setCuRRR <N>                           direct integer RRR (>= 1)
//     /bbr/det/setCuStageT <T> K                      temperature stage [K]
//   from TestWorldMessenger (this example):
//     /bbr/testworld/roundGap <bool>                  also place the round gap (default false)
class TestWorldDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  TestWorldDetectorConstruction();
  ~TestWorldDetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;
  void SetRoundGap(G4bool on) { fRoundGap = on; }

  // The opt-in round gap, after HFSS cylindrical2. The Geant4 radius is the HFSS
  // radius plus 1 um, so the HFSS exit grid's rim points lie strictly inside the
  // solid. Construct() builds the geometry from these and the messenger's
  // guidance prints them.
  static constexpr G4double kRoundGapLength = 0.4 * CLHEP::mm;        // along x, = plate thickness
  static constexpr G4double kRoundGapHfssRadius = 0.050 * CLHEP::mm;
  static constexpr G4double kRoundGapRadius = 0.051 * CLHEP::mm;      // Geant4 solid
  static constexpr G4double kRoundGapPlateZ = -80. * CLHEP::mm;       // plate centre z
  static constexpr G4double kRoundGapPlateHalfWidth = 5. * CLHEP::mm; // half-extent in y and z

 private:
  G4bool fRoundGap = false;
  TestWorldMessenger* fMessenger = nullptr;
};

#endif
