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

  // Two touching Cu plates, each as thick as its crack (OR-020), so every crack volume has exactly
  // the HFSS model's size: crack1 1.0 x 10 x 0.050 mm, crack2 1.5 x 10 x 0.100 mm (propagation,
  // long, gap). Together they keep the old slab's 50 x 50 mm front face at x = 0 (z from -25 to
  // 1.5 mm for plate 1, 1.5 to 25 mm for plate 2). crack1 is centred at (0.5, 0, 0) mm and crack2 at
  // (0.75, 0, 3) mm, as before.
  auto* plate1 = new G4LogicalVolume(new G4Box("solid-CuPlateCrack1", 0.5*mm, 25.*mm, 13.25*mm),
                                     cuMat, "logic-CuPlateCrack1");
  new G4PVPlacement(nullptr, G4ThreeVector(0.5*mm, 0., -11.75*mm),
                    plate1, "CuPlateCrack1", worldLogical, false, 0, true);
  auto* plate2 = new G4LogicalVolume(new G4Box("solid-CuPlateCrack2", 0.75*mm, 25.*mm, 11.75*mm),
                                     cuMat, "logic-CuPlateCrack2");
  new G4PVPlacement(nullptr, G4ThreeVector(0.75*mm, 0., 13.25*mm),
                    plate2, "CuPlateCrack2", worldLogical, false, 0, true);

  // Crack volume names are the HFSS dataset IDs. The data directory adds the
  // frequency: <data root>/waveguides/<id>_<freq>GHz_Ephi={0,1}. BBRCrackLibrary
  // discovers the available frequencies and picks the one nearest each photon.
  // crack1: 50 µm gap (25 µm half-gap), fills plate 1 along x; world z = 0
  {
    const G4String kId = "InfParallelPlate_crack1Rohan";
    auto* solid   = new G4Box(kId, 0.5*mm, 5.*mm, 0.025*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 11.75*mm),
                      logical, kId, plate1, false, 0, true);
  }

  // crack2: 100 µm gap (50 µm half-gap), fills plate 2 along x; world z = 3 mm
  {
    const G4String kId = "InfParallelPlate_crack2";
    auto* solid   = new G4Box(kId, 0.75*mm, 5.*mm, 0.05*mm);
    auto* logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -10.25*mm),
                      logical, kId, plate2, false, 0, true);
  }

  // Straight round gap (HFSS cylindrical2: radius 50 um, length 0.4 mm), only
  // with /bbr/testworld/roundGap true, so the default world and the pinned
  // fixed-seed numbers are unchanged. Its 0.4 mm Cu plate sits at z = -80 mm,
  // clear of the crack plates (|z| <= 25 mm). The vacuum_wg solid is a G4Tubs turned by
  // Ry(+90 deg) inside a G4DisplacedSolid, so its local x is the axis, as
  // BBSimOpBoundaryProcess requires (a bare G4Tubs has it on z: F12 refuses it).
  // The radius is the HFSS radius, 50 um, so the solid equals the HFSS section (OR-020); the
  // exit grid's rim points are clamped 1 nm inside it at load.
  if (fRoundGap) {
    const G4double length = kRoundGapLength, halfWidth = kRoundGapPlateHalfWidth;
    auto* plate = new G4Box("solid-RoundGapPlate", 0.5 * length, halfWidth, halfWidth);
    auto* plateLV = new G4LogicalVolume(plate, cuMat, "logic-RoundGapPlate");
    new G4PVPlacement(nullptr, G4ThreeVector(0.5 * length, 0., kRoundGapPlateZ),
                      plateLV, "RoundGapPlate", worldLogical, false, 0, true);
    const G4String kId = "RoundGap_r50um";
    auto* tube = new G4Tubs("solid-" + kId + "-tube", 0., kRoundGapRadius, 0.5 * length, 0.,
                            CLHEP::twopi);
    G4RotationMatrix ry;
    ry.rotateY(90. * deg);
    auto* gap = new G4DisplacedSolid(kId, tube, G4Transform3D(ry, G4ThreeVector()));
    auto* gapLV = new G4LogicalVolume(gap, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(), gapLV, kId, plateLV, false, 0, true);
  }

  return new G4PVPlacement(nullptr, G4ThreeVector(),
                           worldLogical, "World", nullptr, false, 0, true);
}
