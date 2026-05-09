#ifndef BBRTestDetectorConstruction_hh
#define BBRTestDetectorConstruction_hh
#include "G4VUserDetectorConstruction.hh"

class BBRTestDetectorConstruction : public G4VUserDetectorConstruction {
public:
  BBRTestDetectorConstruction()           = default;
  ~BBRTestDetectorConstruction() override = default;
  G4VPhysicalVolume* Construct() override;
};
#endif
