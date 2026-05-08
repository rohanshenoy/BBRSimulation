#include "BBRPlanckDetectorConstruction.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* BBRPlanckDetectorConstruction::Construct()
{
  auto* solid = new G4Box("World", 10.*cm, 10.*cm, 10.*cm);
  auto* lv    = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), "World");
  return new G4PVPlacement(nullptr, G4ThreeVector(), lv, "World", nullptr, false, 0, true);
}
