#include "TestWorldDetectorConstruction.hh"
#include "TestWorldMessenger.hh"
#include "BBRConfigManager.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"

TestWorldDetectorConstruction::TestWorldDetectorConstruction()
  : fMessenger(new TestWorldMessenger(this)) {}

TestWorldDetectorConstruction::~TestWorldDetectorConstruction() { delete fMessenger; }

G4VPhysicalVolume* TestWorldDetectorConstruction::Construct()
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

  // Crack volume names are the HFSS dataset IDs. The data directory adds the
  // frequency: <data root>/waveguides/<id>_<freq>GHz_Ephi={0,1}. BBRCrackLibrary
  // discovers the available frequencies and picks the one nearest each photon.
  // crack1: 52 µm gap (b=26µm half-width), full-span daughter of CuSlab
  {
    const G4String kId = "InfParallelPlate_crack1Rohan";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.026*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
                      logical, kId, cuLogical, false, 0, true);
  }

  // crack2: 102 µm gap, placed at z=3mm inside CuSlab
  {
    const G4String kId = "InfParallelPlate_crack2";
    auto* solid   = new G4Box(kId, 2.*mm, 5.1*mm, 0.051*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 3.*mm),
                      logical, kId, cuLogical, false, 0, true);
  }

  // Straight round gap (HFSS cylindrical2: radius 50 um, length 0.4 mm), only
  // with /bbr/testworld/roundGap true, so the default world and the pinned
  // fixed-seed numbers are unchanged. Its 0.4 mm Cu plate sits at z = -80 mm,
  // clear of the slab (|z| <= 25 mm). The vacuum_wg solid is a G4Tubs turned by
  // Ry(+90 deg) inside a G4DisplacedSolid, so its local x is the axis, as
  // BBSimOpBoundaryProcess requires (a bare G4Tubs has it on z: F12 refuses it).
  // Radius 51 um = HFSS 50 um + 1 um, so the exit grid's rim points lie strictly
  // inside (the recorded deviation, like the 52/102 um slab gaps).
  if (fRoundGap) {
    const G4double length = 0.4 * mm, radius = 0.051 * mm;
    auto* plate = new G4Box("solid-RoundGapPlate", 0.5 * length, 5. * mm, 5. * mm);
    auto* plateLV = new G4LogicalVolume(plate, cuMat, "logic-RoundGapPlate");
    new G4PVPlacement(nullptr, G4ThreeVector(0.5 * length, 0., -80. * mm),
                      plateLV, "RoundGapPlate", worldLogical, false, 0, true);
    const G4String kId = "RoundGap_r50um";
    auto* tube = new G4Tubs("solid-" + kId + "-tube", 0., radius, 0.5 * length, 0., CLHEP::twopi);
    G4RotationMatrix ry;
    ry.rotateY(90. * deg);
    auto* gap = new G4DisplacedSolid(kId, tube, G4Transform3D(ry, G4ThreeVector()));
    auto* gapLV = new G4LogicalVolume(gap, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(), gapLV, kId, plateLV, false, 0, true);
  }

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLogical, "World", nullptr, false, 0, true);
}
