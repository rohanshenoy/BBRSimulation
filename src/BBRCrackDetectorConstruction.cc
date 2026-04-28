#include "BBRCrackDetectorConstruction.hh"
#include "BBRMaterials.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BBRCrackDetectorConstruction::Construct()
{
  // Vacuum with RINDEX=1 so optical photons propagate freely.
  G4Material* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto mpt = new G4MaterialPropertiesTable();
  G4double energies[] = {1e-6*eV, 1.0*eV};
  G4double rindex[]   = {1.0,     1.0};
  mpt->AddProperty("RINDEX", energies, rindex, 2);
  vac->SetMaterialPropertiesTable(mpt);

  // World: 200 mm cube
  auto worldSolid    = new G4Box("World", 100.*mm, 100.*mm, 100.*mm);
  auto worldLogical  = new G4LogicalVolume(worldSolid, vac, "World");
  auto worldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(),
                                         worldLogical, "World",
                                         nullptr, false, 0, true);

  // Two crack slabs placed side by side (same x-depth, separated in z).
  // Each is named by its HFSS dataset ID; BBSimOpBoundaryProcess auto-loads
  // the matching CSV from HFSSSimData/<name>_Ephi={0,1}/ by volume name.
  // Material "vacuum_wg" marks the volume as a diffraction boundary.
  //
  // Local coordinate convention (identity rotation → world axes):
  //   local +X = propagation (normal_hat)  ← HFSS Z axis
  //   local +Y = long dimension (theta_hat) ← HFSS Y axis
  //   local +Z = gap (phi_hat)              ← HFSS X axis
  //
  // halfY margin (+0.1 mm) and halfZ margin (+1 µm) ensure HFSS exit positions
  // land strictly inside the volume so Geant4 navigation is well-defined.
  //
  // crack1: HFSS ZSize=1 mm depth, XSize=0.05 mm gap, YSize=10 mm
  //         centered at (x=0, y=0, z=0): photon at z=0 hits this crack.
  {
    static const char* kId = "InfParallelPlate_crack1Rohan_500GHz";
    auto solid   = new G4Box(kId, 0.5*mm, 5.1*mm, 0.026*mm);
    auto logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.),
                      logical, kId, worldLogical, false, 0, true);
  }

  // crack2: HFSS ZSize=1.5 mm depth, XSize=0.1 mm gap, YSize=10 mm
  //         centered at (x=0, y=0, z=3 mm): photon offset to z=3 mm hits this crack.
  {
    static const char* kId = "InfParallelPlate_crack2_500GHz";
    auto solid   = new G4Box(kId, 0.751*mm, 5.1*mm, 0.051*mm);
    auto logical = new G4LogicalVolume(solid, BBRMaterials::GetVacuumWG(), kId);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 3.*mm),
                      logical, kId, worldLogical, false, 0, true);
  }

  return worldPhysical;
}
