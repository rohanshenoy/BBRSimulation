#include "BBRTestDetectorConstruction.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

BBRTestDetectorConstruction::BBRTestDetectorConstruction()
  : G4VUserDetectorConstruction()
  , fCuMat(BBRMaterials::GetCopperByName("OFHC_Cu"))
{
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/bbr/det/", "BBR detector geometry commands");
  fMessenger->DeclareMethod("setCuMaterial",
                             &BBRTestDetectorConstruction::SetCuMaterial,
                             "Set Cu wall material: OFHC_Cu | OF_Cu | HP_Cu")
            .SetParameterName("name", true)
            .SetDefaultValue("OFHC_Cu")
            .SetStates(G4State_PreInit);
}

void BBRTestDetectorConstruction::SetCuMaterial(const G4String& name)
{
  G4Material* mat = BBRMaterials::GetCopperByName(name);
  if (!mat) {
    G4cerr << "BBRTestDetectorConstruction::SetCuMaterial: unknown name '"
           << name << "'. Valid: OFHC_Cu, OF_Cu, HP_Cu." << G4endl;
    return;
  }
  fCuMat = mat;
}

G4VPhysicalVolume* BBRTestDetectorConstruction::Construct()
{
  // World: 50 cm cube of G4_Galactic with RINDEX=1
  G4Material* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  {
    const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
    const std::vector<G4double> ri = {1., 1.};
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, ri);
    vac->SetMaterialPropertiesTable(mpt);
  }

  auto* worldSolid   = new G4Box("solid-World",  250.*mm, 250.*mm, 250.*mm);
  auto* worldLogical = new G4LogicalVolume(worldSolid, vac, "logic-World");

  // Cu wall: center at (2mm, 0, 0), front face at x=0, back face at x=4mm
  auto* cuSolid   = new G4Box("solid-CuSlab", 2.*mm, 25.*mm, 25.*mm);
  auto* cuLogical = new G4LogicalVolume(cuSolid, fCuMat, "logic-CuSlab");
  new G4PVPlacement(nullptr, G4ThreeVector(2.*mm, 0., 0.),
                    cuLogical, "CuSlab", worldLogical, false, 0, true);

  // crack1: 52 µm gap (b=26µm half-width), full-span daughter of CuSlab
  {
    const G4String kId = "InfParallelPlate_crack1Rohan_500GHz";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.026*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
                      logical, kId, cuLogical, false, 0, true);
  }

  // crack2: 102 µm gap, placed at z=3mm inside CuSlab
  {
    const G4String kId = "InfParallelPlate_crack2_500GHz";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.051*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 3.*mm),
                      logical, kId, cuLogical, false, 0, true);
  }

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLogical, "World", nullptr, false, 0, true);
}
