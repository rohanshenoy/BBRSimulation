#include "BBRReflectanceDetectorConstruction.hh"
#include "BBRMaterials.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BBRReflectanceDetectorConstruction::Construct()
{
  auto* worldBox = new G4Box("World", 50.*mm, 50.*mm, 50.*mm);
  auto* worldLV  = new G4LogicalVolume(worldBox, BBRMaterials::GetVacuumWG(), "World");

  // Cu slab: front (−x) face at x=0, centre at x=+1 mm.
  auto* cuBox = new G4Box("CuSlab", 1.*mm, 5.*mm, 5.*mm);
  auto* cuLV  = new G4LogicalVolume(cuBox, BBRMaterials::GetOFHCCopper(), "CuSlab");
  new G4PVPlacement(nullptr, G4ThreeVector(1.*mm, 0., 0.),
                    cuLV, "CuSlab", worldLV, false, 0, true);

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLV, "World", nullptr, false, 0, true);
}
