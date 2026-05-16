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
