#ifndef BBRMaterials_hh
#define BBRMaterials_hh

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>
#include <vector>

namespace BBRMaterials {

// ---------------------------------------------------------------------------
// Internal helper — not part of the public API.
// Builds a single-element material with Hagen-Rubens reflectance table.
//   name             : G4Material name (checked for existence first)
//   Z, A_g_mol       : atomic number and mass [g/mol]
//   density_g_cm3    : density [g/cm³]
//   sigma_SI         : electrical conductivity at 4 K [S/m]
// ---------------------------------------------------------------------------
inline G4Material* BuildHagRubMaterial(const G4String& name,
                                        G4double Z,
                                        G4double A_g_mol,
                                        G4double density_g_cm3,
                                        G4double sigma_SI)
{
  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat) return mat;
  mat = new G4Material(name, Z, A_g_mol*g/mole, density_g_cm3*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;         // 50 GHz
  const G4double Emax  = 8.27e-2*eV;         // 20 THz
  const G4double eps0  = 8.8541878128e-12;   // F/m (SI)
  const G4double h_eVs = 4.13566769692e-15;  // eV·s (Planck constant)

  std::vector<G4double> energies(N), refls(N);
  const G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma_SI);
    refls[i]       = std::max(0., std::min(1., R));
  }
  const std::vector<G4double> e2 = {Emin, Emax};
  const std::vector<G4double> ri = {1., 1.};

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,      ri);
  mpt->AddProperty("REFLECTIVITY", energies, refls);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

// Near-vacuum material that flags crack volumes.
// RINDEX=1, no absorption, distinct name detected by BBSimOpBoundaryProcess.
inline G4Material* GetVacuumWG()
{
  G4Material* mat = G4Material::GetMaterial("vacuum_wg", false);
  if (mat) return mat;
  mat = new G4Material("vacuum_wg", 1., 1.008*g/mole, 1e-25*g/cm3);
  const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
  const std::vector<G4double> ri = {1.0,     1.0};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", e, ri);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// R=0 everywhere: every photon absorbed on contact.
inline G4Material* GetPerfectAbsorber()
{
  G4Material* mat = G4Material::GetMaterial("BBR_PerfectAbsorber", false);
  if (mat) return mat;
  mat = new G4Material("BBR_PerfectAbsorber", 1., 1.008*g/mole, 1e-25*g/cm3);
  const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
  const std::vector<G4double> ri = {1., 1.};
  const std::vector<G4double> r  = {0., 0.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e, ri);
  mpt->AddProperty("REFLECTIVITY", e, r);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// R=1 everywhere: every photon reflects specularly.
inline G4Material* GetPerfectReflector()
{
  G4Material* mat = G4Material::GetMaterial("BBR_PerfectReflector", false);
  if (mat) return mat;
  mat = new G4Material("BBR_PerfectReflector", 1., 1.008*g/mole, 1e-25*g/cm3);
  const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
  const std::vector<G4double> ri = {1., 1.};
  const std::vector<G4double> r  = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e, ri);
  mpt->AddProperty("REFLECTIVITY", e, r);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// OFHC Cu at 4 K, RRR=100. σ = 100 × σ_RT(Cu) = 5.96×10⁹ S/m.
inline G4Material* GetOFHCCopper()
{ return BuildHagRubMaterial("OFHC_Cu", 29., 63.546, 8.96, 5.96e9); }

// OF copper at 4 K (Serov 2016, Fig. 6). σ_eff = 1/ρ₀, ρ₀=0.56×10⁻⁸ Ω·m.
inline G4Material* GetOFCopperSerov()
{ return BuildHagRubMaterial("OF_Cu",   29., 63.546, 8.96, 1.786e8); }

// HP (hydrogen-annealed) copper at 4 K (Serov 2016, Fig. 8).
// σ_eff back-calculated from D=0.55×10⁻³ at 230 GHz.
inline G4Material* GetHPCopperSerov()
{ return BuildHagRubMaterial("HP_Cu",   29., 63.546, 8.96, 3.383e8); }

// Lookup by name string — used by BBRTestDetectorConstruction messenger.
inline G4Material* GetCopperByName(const G4String& name)
{
  if (name == "OFHC_Cu") return GetOFHCCopper();
  if (name == "OF_Cu")   return GetOFCopperSerov();
  if (name == "HP_Cu")   return GetHPCopperSerov();
  return nullptr;
}

} // namespace BBRMaterials

#endif
