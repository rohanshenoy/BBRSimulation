#include "BBRTestDetectorConstruction.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BBRTestDetectorConstruction::Construct()
{
  // World: 50 cm cube of G4_Galactic with RINDEX=1
  G4Material* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto* worldMPT = new G4MaterialPropertiesTable();
  G4double e[]  = {1e-6*eV, 1.0*eV};
  G4double ri[] = {1., 1.};
  worldMPT->AddProperty("RINDEX", e, ri, 2);
  vac->SetMaterialPropertiesTable(worldMPT);

  auto* worldBox = new G4Box("World", 250.*mm, 250.*mm, 250.*mm);
  auto* worldLV  = new G4LogicalVolume(worldBox, vac, "World");

  // Cu wall: center at (2mm,0,0), front face at x=0, back face at x=4mm
  auto* cuBox = new G4Box("CuSlab", 2.*mm, 25.*mm, 25.*mm);
  auto* cuLV  = new G4LogicalVolume(cuBox, BBRMaterials::GetOFHCCopper(), "CuSlab");
  new G4PVPlacement(nullptr, G4ThreeVector(2.*mm, 0., 0.),
                    cuLV, "CuSlab", worldLV, false, 0, true);

  // crack1: daughter of CuSlab in Cu local frame
  {
    static const char* kId = "InfParallelPlate_crack1Rohan_500GHz";
    auto* solid   = new G4Box(kId, 0.5*mm, 5.1*mm, 0.026*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(-1.499*mm, 0., 0.),
                      logical, kId, cuLV, false, 0, true);
  }

  // crack2: daughter of CuSlab in Cu local frame
  {
    static const char* kId = "InfParallelPlate_crack2_500GHz";
    auto* solid   = new G4Box(kId, 0.751*mm, 5.1*mm, 0.051*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(-1.248*mm, 0., 3.*mm),
                      logical, kId, cuLV, false, 0, true);
  }

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLV, "World", nullptr, false, 0, true);
}
