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
#include "G4VSolid.hh"
#include "CADMesh.hh"
#include <cmath>
#include <fstream>

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

  // Lay the tube axis (local z) along world +x with the warm aperture at
  // x = -50 mm. The tube is symmetric under z -> -z, so the rotation sign is
  // immaterial.
  auto* rot = new G4RotationMatrix();
  rot->rotateY(90.*deg);
  const G4double xTubeMin = -50.*mm;
  const G4double xCenter  = xTubeMin + fLength / 2.;

  new G4PVPlacement(rot, G4ThreeVector(xCenter, 0., 0.), tubeLV,
                    "LightPipeWall", worldLV, false, 0, true);

  // BBRLightPipe reuses BBRTestPGA, whose Planck emitter box is set with
  // /bbr/thermal/emitterCenter and /bbr/thermal/emitterSize. The test-world
  // default (1x20x20 mm centred at x = -50 mm) overlaps this wall: its
  // emitting face sits 0.5 mm inside the tube and reaches r = 14 mm, so
  // primaries would be created inside the copper and pass through it. Warn
  // when the configured emitter reaches into the wall region.
  const G4ThreeVector ec = BBRConfigManager::GetEmitterCenter_mm() * mm;
  const G4ThreeVector es = BBRConfigManager::GetEmitterSize_mm() * mm;
  const G4double emitterXMax = ec.x() + 0.5 * es.x();
  const G4double emitterRMax = std::hypot(ec.y(), ec.z())
                             + std::hypot(0.5 * es.y(), 0.5 * es.z());
  if (emitterXMax > xTubeMin && emitterRMax > fBore) {
    G4ExceptionDescription ed;
    ed << "Planck emitter box (centre " << ec / mm << " mm, extents " << es / mm
       << " mm) reaches into the light-pipe wall (x >= " << xTubeMin / mm
       << " mm, r in [" << fBore / mm << ", " << (fBore + fWallThickness) / mm
       << "] mm): primaries would be created inside copper. Set "
          "/bbr/thermal/emitterCenter and /bbr/thermal/emitterSize so the "
          "emitter sits upstream of the aperture and inside the bore "
          "(see lightpipe.mac).";
    G4Exception("BBRLightPipeDetectorConstruction::BuildParametric", "LP002",
                JustWarning, ed);
  }

  G4cout << "[BBR] LightPipe parametric: bore=" << fBore/mm
         << "mm length=" << fLength/mm << "mm wall=" << fWallThickness/mm
         << "mm material=" << wall->GetName() << G4endl;
}

void BBRLightPipeDetectorConstruction::BuildFromCAD(G4LogicalVolume* worldLV)
{
  if (fStlPath.empty()) {
    G4Exception("BBRLightPipeDetectorConstruction::BuildFromCAD", "LP010",
                FatalException, "cad mode requires /bbr/lightpipe/stlPath");
  }
  // Fail loudly on a missing file rather than building an empty world.
  if (!std::ifstream(fStlPath).good()) {
    G4Exception("BBRLightPipeDetectorConstruction::BuildFromCAD", "LP011",
                FatalException, ("STL not found: " + fStlPath).c_str());
  }

  G4VSolid* solid = nullptr;
  try {
    auto mesh = CADMesh::TessellatedMesh::FromSTL(fStlPath);
    mesh->SetScale(mm);
    solid = mesh->GetSolid();
  } catch (...) {
    G4Exception("BBRLightPipeDetectorConstruction::BuildFromCAD", "LP012",
                FatalException,
                ("failed to parse STL (built-in reader is ASCII-only): "
                 + fStlPath).c_str());
  }

  auto* lv = new G4LogicalVolume(solid, ResolveWallMaterial(), "logic-LightPipe");
  new G4PVPlacement(nullptr, G4ThreeVector(), lv, "LightPipeWall",
                    worldLV, false, 0, true);

  G4cout << "[BBR] LightPipe CAD: " << fStlPath
         << " material=" << lv->GetMaterial()->GetName() << G4endl;
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

  if (fMode == "cad") BuildFromCAD(worldLogical);
  else                BuildParametric(worldLogical);
  return worldPhys;
}
