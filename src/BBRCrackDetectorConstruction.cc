#include "BBRCrackDetectorConstruction.hh"

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

  // TEM_waveguide slab at origin.
  // HFSS model: propagation along +x (normal_hat), gap b=0.05 mm along +z (phi_hat),
  // long dimension along +y (theta_hat).
  //   halfX = a/2 = 5 mm   (crack depth)
  //   halfY = 50 mm         (long dimension, truncated from infinite)
  //   halfZ = b/2 = 0.025 mm (gap half-width)
  auto crackSolid   = new G4Box("TEM_waveguide_crack", 5.*mm, 50.*mm, 0.025*mm);
  auto crackLogical = new G4LogicalVolume(crackSolid, vac, "TEM_waveguide_crack");
  new G4PVPlacement(nullptr, G4ThreeVector(),
                    crackLogical, "TEM_waveguide_crack",
                    worldLogical, false, 0, true);

  return worldPhysical;
}
