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

  // Waveguide crack slab at origin. Material "vacuum_wg" marks it as a crack;
  // BBSimOpBoundaryProcess detects the vacuum→vacuum_wg boundary and routes
  // photons through the HFSS diffraction model.
  //
  // Local coordinate convention: local +X = propagation (normal_hat),
  //   local +Y = long dimension (theta_hat), local +Z = gap (phi_hat).
  // With identity rotation: local axes align with world +x/+y/+z.
  //   halfX = 5 mm    (crack depth along propagation axis)
  //   halfY = 50 mm   (long dimension a, truncated from infinite)
  //   halfZ = 0.025 mm (gap half-width b/2 = 25 µm)
  //
  // Volume name = HFSS dataset ID (path prefix under HFSSSimData/).
  static const char* kDatasetId = "InfParallelPlate_crack1Rohan_500GHz";
  auto crackSolid   = new G4Box(kDatasetId, 5.*mm, 50.*mm, 0.025*mm);
  auto crackLogical = new G4LogicalVolume(crackSolid, BBRMaterials::GetVacuumWG(), kDatasetId);
  new G4PVPlacement(nullptr, G4ThreeVector(),
                    crackLogical, kDatasetId,
                    worldLogical, false, 0, true);

  return worldPhysical;
}
