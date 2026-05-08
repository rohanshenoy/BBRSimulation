#ifndef BBRPlanckDetectorConstruction_hh
#define BBRPlanckDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

// Empty 20cm vacuum_wg cube. Photons propagate freely until they exit the world.
class BBRPlanckDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override;
};

#endif
