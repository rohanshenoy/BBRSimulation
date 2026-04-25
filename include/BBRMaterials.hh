#ifndef BBRMaterials_hh
#define BBRMaterials_hh

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"

// Seed of the BBR material database. Provides shared material definitions
// used across BBR detector constructions.
namespace BBRMaterials {

// Returns the "vacuum_wg" material: optically identical to G4_Galactic (RINDEX=1,
// no absorption) but with a distinct name. Any volume filled with this material
// is treated as a waveguide crack by BBSimOpBoundaryProcess.
// Singleton: created on first call, reused thereafter.
inline G4Material* GetVacuumWG()
{
  G4Material* mat = G4Material::GetMaterial("vacuum_wg", false);
  if (mat) return mat;
  mat = new G4Material("vacuum_wg", 1., 1.008 * g / mole, 1e-25 * g / cm3);
  auto* mpt = new G4MaterialPropertiesTable();
  G4double energies[] = {1e-6 * eV, 1.0 * eV};
  G4double rindex[]   = {1.0,       1.0};
  mpt->AddProperty("RINDEX", energies, rindex, 2);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

} // namespace BBRMaterials

#endif
