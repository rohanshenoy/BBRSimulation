#include "BBRReflectanceDetectorConstruction.hh"
#include "BBRMaterials.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BBRReflectanceDetectorConstruction::Construct()
{
  // Free-space world: G4_Galactic with RINDEX=1, same pattern as
  // BBRCrackDetectorConstruction.  vacuum_wg is reserved for crack volumes
  // only; using it here would confuse the diffraction-boundary dispatch.
  G4Material* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto* worldMPT = new G4MaterialPropertiesTable();
  G4double e[]  = {1e-6*eV, 1.0*eV};
  G4double ri[] = {1., 1.};
  worldMPT->AddProperty("RINDEX", e, ri, 2);
  vac->SetMaterialPropertiesTable(worldMPT);

  auto* worldBox = new G4Box("World", 50.*mm, 50.*mm, 50.*mm);
  auto* worldLV  = new G4LogicalVolume(worldBox, vac, "World");

  // Cu slab: front (−x) face at x=0, centre at x=+1 mm.
  auto* cuBox = new G4Box("CuSlab", 1.*mm, 5.*mm, 5.*mm);
  auto* cuLV  = new G4LogicalVolume(cuBox, BBRMaterials::GetOFHCCopper(), "CuSlab");
  new G4PVPlacement(nullptr, G4ThreeVector(1.*mm, 0., 0.),
                    cuLV, "CuSlab", worldLV, false, 0, true);

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLV, "World", nullptr, false, 0, true);
}
