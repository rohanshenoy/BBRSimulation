#ifndef BBRTestDetectorConstruction_hh
#define BBRTestDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4Material.hh"
#include <memory>

// Test geometry: 50cm world, 4mm Cu slab at x=2mm with two vacuum_wg crack daughters.
// Default Cu material: OFHC_Cu (RRR=100 Hagen-Rubens).
// Change via: /bbr/det/setCuMaterial <OFHC_Cu|OF_Cu|HP_Cu>  (before /run/initialize)
class BBRTestDetectorConstruction : public G4VUserDetectorConstruction {
public:
  BBRTestDetectorConstruction();
  ~BBRTestDetectorConstruction() override = default;
  G4VPhysicalVolume* Construct() override;

  // Called by messenger; also callable directly from main for testing.
  void SetCuMaterial(const G4String& name);

private:
  G4Material*                          fCuMat;
  std::unique_ptr<G4GenericMessenger>  fMessenger;
};

#endif
