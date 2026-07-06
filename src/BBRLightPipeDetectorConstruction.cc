#include "BBRLightPipeDetectorConstruction.hh"
#include "BBRLightPipeMessenger.hh"
#include "BBRConfigManager.hh"
#include "BBRMaterials.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

BBRLightPipeDetectorConstruction::BBRLightPipeDetectorConstruction()
  : fMessenger(new BBRLightPipeMessenger(this)) {}

BBRLightPipeDetectorConstruction::~BBRLightPipeDetectorConstruction()
{
  delete fMessenger;
}

G4Material* BBRLightPipeDetectorConstruction::ResolveWallMaterial()
{
  if (fWallMaterialName == "reflector")
    return BBRMaterials::GetPerfectReflector();
  // Default: Cu with RRR / stage-T from the shared run config.
  return BBRMaterials::GetCopper(BBRConfigManager::GetCuRRR(),
                                 BBRConfigManager::GetCuStageT_K());
}

void BBRLightPipeDetectorConstruction::BuildParametric(G4LogicalVolume* worldLV)
{
  if (fBore <= 0. || fLength <= 0. || fWallThickness <= 0.) {
    G4Exception("BBRLightPipeDetectorConstruction::BuildParametric", "LP001",
                FatalException, "bore, length, wallThickness must all be > 0");
  }

  G4Material* wall = ResolveWallMaterial();

  // Hollow tube: bore is the mother (vacuum) volume; only the wall is a solid.
  auto* tube = new G4Tubs("solid-LightPipe", fBore, fBore + fWallThickness,
                          fLength / 2., 0., CLHEP::twopi);
  auto* tubeLV = new G4LogicalVolume(tube, wall, "logic-LightPipe");

  // Lay the tube axis (local z) along world +x; warm aperture just past the
  // emitter patch at x = -50 mm so BBRTestPGA illuminates the bore. The tube is
  // symmetric under z -> -z, so the rotation sign is immaterial.
  auto* rot = new G4RotationMatrix();
  rot->rotateY(90.*deg);
  const G4double xCenter = -50.*mm + fLength / 2.;

  new G4PVPlacement(rot, G4ThreeVector(xCenter, 0., 0.), tubeLV,
                    "LightPipeWall", worldLV, false, 0, true);

  G4cout << "[BBR] LightPipe parametric: bore=" << fBore/mm
         << "mm length=" << fLength/mm << "mm wall=" << fWallThickness/mm
         << "mm material=" << wall->GetName() << G4endl;
}

G4VPhysicalVolume* BBRLightPipeDetectorConstruction::Construct()
{
  // Vacuum world (matches the test world: G4_Galactic, RINDEX=1).
  G4Material* vac =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  if (!vac->GetMaterialPropertiesTable()) {
    const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
    const std::vector<G4double> ri = {1., 1.};
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, ri);
    vac->SetMaterialPropertiesTable(mpt);
  }

  auto* worldSolid   = new G4Box("solid-World", 250.*mm, 250.*mm, 250.*mm);
  auto* worldLogical = new G4LogicalVolume(worldSolid, vac, "logic-World");
  auto* worldPhys = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical,
                                      "World", nullptr, false, 0, true);

  BuildParametric(worldLogical);
  return worldPhys;
}
