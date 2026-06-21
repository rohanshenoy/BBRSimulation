#include "BBRTestDetectorConstruction.hh"
#include "BBRConfigManager.hh"
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
  const G4int    rrr     = BBRConfigManager::GetCuRRR();
  const G4double stageT  = BBRConfigManager::GetCuStageT_K();
  G4Material* cuMat = BBRMaterials::GetCopper(rrr, stageT);
  G4cout << "[BBR] Cu wall material: " << cuMat->GetName()
         << "  (RRR=" << rrr << ", T=" << stageT << " K)" << G4endl;

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
  auto* cuLogical = new G4LogicalVolume(cuSolid, cuMat, "logic-CuSlab");
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
