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
{
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/bbr/det/", "BBR detector geometry commands");

  fMessenger->DeclareMethod(
      "setCuMaterial",
      &BBRTestDetectorConstruction::SetCuMaterial,
      "Named Cu alias: OFHC_Cu (RRR=100) | OF_Cu (RRR=3) | HP_Cu (RRR=6). "
      "Sets RRR and resets temperature to 4 K.")
    .SetParameterName("alias", true)
    .SetDefaultValue("OFHC_Cu")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);  // geometry: built on master, not broadcast to workers

  fMessenger->DeclareMethod(
      "setCuRRR",
      &BBRTestDetectorConstruction::SetCuRRR,
      "Cu Residual Resistance Ratio (integer >= 1). "
      "Derives sigma_DC = RRR * 5.96e7 S/m. "
      "Combine with setCuStageT to set the temperature stage. "
      "Must be issued before /run/initialize.")
    .SetParameterName("RRR", false)
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);  // geometry: built on master, not broadcast to workers

  fMessenger->DeclareMethodWithUnit(
      "setCuStageT", "kelvin",
      &BBRTestDetectorConstruction::SetCuStageT,
      "Temperature stage [K] for the Cu reflectance table (default 4 K). "
      "Below 50 K: sigma = RRR * sigma_RT. "
      "Above 50 K: Matthiessen's rule (1/T phonon approximation). "
      "Must be issued before /run/initialize.")
    .SetParameterName("T", false)
    .SetRange("T > 0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);  // geometry: built on master, not broadcast to workers
}

void BBRTestDetectorConstruction::SetCuMaterial(const G4String& alias)
{
  if      (alias == "OFHC_Cu") { fCuRRR = 100; fStageT_K = 4.0; }
  else if (alias == "OF_Cu")   { fCuRRR =   3; fStageT_K = 4.0; }
  else if (alias == "HP_Cu")   { fCuRRR =   6; fStageT_K = 4.0; }
  else {
    G4cerr << "[BBR] setCuMaterial: unknown alias '" << alias
           << "'.  Valid: OFHC_Cu, OF_Cu, HP_Cu." << G4endl;
  }
}

void BBRTestDetectorConstruction::SetCuRRR(G4int rrr)
{
  if (rrr < 1) {
    G4cerr << "[BBR] setCuRRR: RRR must be >= 1, got " << rrr << G4endl;
    return;
  }
  fCuRRR = rrr;
}

void BBRTestDetectorConstruction::SetCuStageT(G4double T_K)
{
  if (T_K <= 0.) {
    G4cerr << "[BBR] setCuStageT: temperature must be > 0 K, got " << T_K << G4endl;
    return;
  }
  fStageT_K = T_K;
}

G4VPhysicalVolume* BBRTestDetectorConstruction::Construct()
{
  G4Material* cuMat = BBRMaterials::GetCopper(fCuRRR, fStageT_K);
  G4cout << "[BBR] Cu wall material: " << cuMat->GetName()
         << "  (RRR=" << fCuRRR << ", T=" << fStageT_K << " K)" << G4endl;

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
