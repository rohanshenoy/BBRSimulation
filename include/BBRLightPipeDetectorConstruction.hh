#ifndef BBRLightPipeDetectorConstruction_hh
#define BBRLightPipeDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class BBRLightPipeMessenger;

// 4 K -> mixing-chamber light pipe. Parametric G4Tubs now; CAD .STL in a later
// task. Physics rides on the wall material's REFLECTIVITY MPT, so
// BBSimOpBoundaryProcess handles the wall identically for tube or mesh.
class BBRLightPipeDetectorConstruction : public G4VUserDetectorConstruction {
public:
  BBRLightPipeDetectorConstruction();
  ~BBRLightPipeDetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  // Geometry parameters (set by BBRLightPipeMessenger, read in Construct()).
  void SetBore(G4double v)                { fBore = v; }
  void SetLength(G4double v)              { fLength = v; }
  void SetWallThickness(G4double v)       { fWallThickness = v; }
  void SetWallMaterial(const G4String& m) { fWallMaterialName = m; }

private:
  void        BuildParametric(G4LogicalVolume* worldLV);
  G4Material* ResolveWallMaterial();

  G4double fBore            = 5.*mm;
  G4double fLength          = 100.*mm;
  G4double fWallThickness   = 2.*mm;
  G4String fWallMaterialName = "Cu";

  BBRLightPipeMessenger* fMessenger = nullptr;
};

#endif
