#ifndef BBRMaterials_hh
#define BBRMaterials_hh

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>

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

// OFHC Cu at 4 K, RRR=100.  Hagen-Rubens: R(ν) = 1 − 2·sqrt(2ε₀ω/σ)
// σ = RRR·σ_RT(Cu) = 100 × 5.96e7 S/m = 5.96e9 S/m.
// 20 log-spaced energies, 50 GHz – 20 THz (2.07e-4 – 8.27e-2 eV).
// "REFLECTIVITY" matches YYC's property name; safe on Material2 MPT because
// stock G4OpBoundaryProcess only reads it from G4OpticalSurface, not Material2.
inline G4Material* GetOFHCCopper()
{
  G4Material* mat = G4Material::GetMaterial("OFHC_Cu", false);
  if (mat) return mat;
  mat = new G4Material("OFHC_Cu", 29., 63.546*g/mole, 8.96*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;        // 50 GHz
  const G4double Emax  = 8.27e-2*eV;        // 20 THz
  const G4double eps0  = 8.8541878128e-12;  // F/m
  const G4double sigma = 5.96e9;            // S/m
  const G4double h_eVs = 4.13566769692e-15; // eV·s

  G4double energies[N], refls[N];
  G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma);
    refls[i]       = std::max(0., std::min(1., R));
  }
  G4double e2[]     = {Emin, Emax};
  G4double rindex[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       rindex, 2);
  mpt->AddProperty("REFLECTIVITY", energies, refls,  N);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// R=0 everywhere: every photon is absorbed.
inline G4Material* GetPerfectAbsorber()
{
  G4Material* mat = G4Material::GetMaterial("BBR_PerfectAbsorber", false);
  if (mat) return mat;
  mat = new G4Material("BBR_PerfectAbsorber", 1., 1.008*g/mole, 1e-25*g/cm3);
  G4double e[]  = {1e-6*eV, 1.0*eV};
  G4double r[]  = {0., 0.};
  G4double ri[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e, ri, 2);
  mpt->AddProperty("REFLECTIVITY", e, r,  2);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// R=1 everywhere: every photon reflects specularly.
inline G4Material* GetPerfectReflector()
{
  G4Material* mat = G4Material::GetMaterial("BBR_PerfectReflector", false);
  if (mat) return mat;
  mat = new G4Material("BBR_PerfectReflector", 1., 1.008*g/mole, 1e-25*g/cm3);
  G4double e[]  = {1e-6*eV, 1.0*eV};
  G4double r[]  = {1., 1.};
  G4double ri[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e, ri, 2);
  mpt->AddProperty("REFLECTIVITY", e, r,  2);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// OF copper at 4 K (Serov et al., IEEE TMT 2016, Fig. 6).
// 99.97% Cu, residual-resistance dominated: σ_eff = 1/ρ₀, ρ₀=0.56e-8 Ω·m.
// Hagen-Rubens D(f) ≈ 0.61e-3 at 150 GHz (Serov: 0.58e-3, within 6%).
inline G4Material* GetOFCopperSerov()
{
  G4Material* mat = G4Material::GetMaterial("OF_Cu", false);
  if (mat) return mat;
  mat = new G4Material("OF_Cu", 29., 63.546*g/mole, 8.96*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;        // 50 GHz
  const G4double Emax  = 8.27e-2*eV;        // 20 THz
  const G4double eps0  = 8.8541878128e-12;  // F/m
  const G4double sigma = 1.786e8;           // S/m  (= 1/ρ₀, ρ₀=0.56e-8 Ω·m)
  const G4double h_eVs = 4.13566769692e-15; // eV·s

  G4double energies[N], refls[N];
  G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma);
    refls[i]       = std::max(0., std::min(1., R));
  }
  G4double e2[]     = {Emin, Emax};
  G4double rindex[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       rindex, 2);
  mpt->AddProperty("REFLECTIVITY", energies, refls,  N);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// HP copper (hydrogen-annealed) at 4 K (Serov et al., IEEE TMT 2016, Fig. 8).
// 99.999% Cu, annealed; σ_eff back-calculated from D=0.55e-3 at 230 GHz.
// Hagen-Rubens D(f) ≈ 0.55e-3 at 230 GHz by construction.
inline G4Material* GetHPCopperSerov()
{
  G4Material* mat = G4Material::GetMaterial("HP_Cu", false);
  if (mat) return mat;
  mat = new G4Material("HP_Cu", 29., 63.546*g/mole, 8.96*g/cm3);

  const int      N     = 20;
  const G4double Emin  = 2.07e-4*eV;        // 50 GHz
  const G4double Emax  = 8.27e-2*eV;        // 20 THz
  const G4double eps0  = 8.8541878128e-12;  // F/m
  const G4double sigma = 3.383e8;           // S/m  (back-calc, D=0.55e-3 @ 230 GHz)
  const G4double h_eVs = 4.13566769692e-15; // eV·s

  G4double energies[N], refls[N];
  G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]    = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    G4double nu    = (energies[i]/eV) / h_eVs;
    G4double omega = 2.*CLHEP::pi*nu;
    G4double R     = 1. - 2.*std::sqrt(2.*eps0*omega/sigma);
    refls[i]       = std::max(0., std::min(1., R));
  }
  G4double e2[]     = {Emin, Emax};
  G4double rindex[] = {1., 1.};
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       rindex, 2);
  mpt->AddProperty("REFLECTIVITY", energies, refls,  N);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

} // namespace BBRMaterials

#endif
