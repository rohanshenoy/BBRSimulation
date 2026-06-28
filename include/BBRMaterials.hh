#ifndef BBRMaterials_hh
#define BBRMaterials_hh

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

namespace BBRMaterials {

// ---------------------------------------------------------------------------
// Internal helper — not part of the public API.
// Builds a single-element material with a Drude-model REFLECTIVITY table.
//   name          : G4Material name (returned as-is if it already exists)
//   Z, A_g_mol    : atomic number and mass [g/mol]
//   density_g_cm3 : density [g/cm³]
//   RRR           : Residual Resistance Ratio (material quality parameter)
//   T_K           : temperature [K]
//
// Physics (Griffiths §9.4 generalized to complex σ + Matthiessen's rule):
//   σ_DC = RRR × σ_RT               (T < 50 K, phonons frozen out)
//   τ    = σ_DC × mₑ / (nₑ e²)
//   σ(ω) = σ_DC / (1 − iωτ)        (full complex Drude AC conductivity)
//   ε̃(ω) = 1 + iσ(ω)/(ε₀ω)         (k̃² = (ω/c)² ε̃;  ñ = √ε̃)
//   R(ω) = |(ñ−1)/(ñ+1)|²          (Fresnel, normal incidence)
//
// Im σ must NOT be dropped: for ωτ ≳ 1 it supplies the −ωp²τ²/(1+ω²τ²)
// plasma term in Re ε̃, which keeps R near 1 (relaxation regime,
// D ≈ 2/(ωp τ)). Using only Re σ in the real-σ Griffiths closed form
// overestimates absorptance by ~30× at 500 GHz (RRR=100, 4 K) and by
// orders of magnitude at 20 THz.
// ---------------------------------------------------------------------------
inline G4Material* BuildDrudeMaterial(const G4String& name,
                                       G4double Z,
                                       G4double A_g_mol,
                                       G4double density_g_cm3,
                                       G4int    RRR,
                                       G4double T_K)
{
  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat) return mat;
  mat = new G4Material(name, Z, A_g_mol*g/mole, density_g_cm3*g/cm3);

  // SI constants — CLHEP values are not SI, so use literals.
  const G4double sigma_RT = 5.96e7;            // S/m, universal for Cu at 273 K
  const G4double n_e      = 8.49e28;           // m^-3, free electron density
  const G4double m_e_kg   = 9.109e-31;         // kg
  const G4double e_C      = 1.602e-19;         // C
  const G4double eps0_SI  = 8.8541878128e-12;  // F/m
  const G4double c_SI     = 2.998e8;           // m/s
  const G4double h_eVs    = 4.13566769692e-15; // eV·s (Planck constant)

  // DC conductivity via Matthiessen's rule.
  // The linear σ_phonon(T) ≈ σ_RT × 273/T is only valid above ~50 K.
  // Below 50 K phonons are frozen out and σ_DC ≈ σ_imp = RRR × σ_RT.
  const G4double sigma_imp = static_cast<G4double>(RRR) * sigma_RT;
  G4double sigma_DC;
  if (T_K >= 50.) {
    const G4double sigma_ph = sigma_RT * 273. / T_K;
    sigma_DC = 1. / (1./sigma_imp + 1./sigma_ph);
  } else {
    sigma_DC = sigma_imp;
  }

  // Drude scattering time: τ = σ_DC mₑ / (nₑ e²)
  const G4double tau = sigma_DC * m_e_kg / (n_e * e_C * e_C);

  // Build 24 log-spaced reflectivity points from 10 GHz to 20 THz.
  // Lower bound matches the Planck-emitter CDF (10 GHz) so sub-50-GHz
  // photons are no longer clamped to the table edge.
  const int      N     = 24;
  const G4double Emin  = 4.14e-5*eV;   // 10 GHz
  const G4double Emax  = 8.27e-2*eV;   // 20 THz

  std::vector<G4double> energies(N), refls(N);
  const G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]          = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    const G4double nu    = (energies[i]/eV) / h_eVs;    // Hz
    const G4double omega = 2. * CLHEP::pi * nu;          // rad/s

    // Full complex Drude conductivity: σ(ω) = σ_DC/(1−iωτ)
    const std::complex<G4double> sigma =
        sigma_DC / std::complex<G4double>(1., -omega*tau);

    // ε̃ = 1 + iσ/(ε₀ω);  ñ = √ε̃
    const std::complex<G4double> eps_t =
        1. + std::complex<G4double>(0., 1.) * sigma / (eps0_SI * omega);
    const std::complex<G4double> n_t = std::sqrt(eps_t);

    // Normal-incidence Fresnel: R = |(ñ−1)/(ñ+1)|²
    const G4double R = std::norm((n_t - 1.) / (n_t + 1.));
    refls[i] = std::max(0., std::min(1., R));
  }

  // RINDEX flat at 1 (two boundary points suffice); REFLECTIVITY uses all N.
  const std::vector<G4double> e2 = {Emin, Emax};
  const std::vector<G4double> ri = {1., 1.};

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       e2,       ri);
  mpt->AddProperty("REFLECTIVITY", energies, refls);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// ---------------------------------------------------------------------------
// Internal helper — not part of the public API.
// Builds a NIST-based dielectric with a flat RINDEX and a loss-tangent-derived
// ABSLENGTH table (bulk absorption handled by stock G4OpAbsorption; the
// vacuum→dielectric Fresnel boundary by stock G4OpBoundaryProcess once RINDEX
// is set — no wrapper interception needed).
//   name      : BBR material name (cached/returned if it already exists)
//   nistBase  : NIST base material to clone (e.g. "G4_KAPTON", "G4_Si")
//   n_index   : real refractive index, flat over the band
//   tan_delta : loss tangent (frequency-independent here — the cryogenic
//               single-fit regime)
//
// Physics: a dielectric with loss tangent tan δ attenuates intensity as
//   α(ν) = 2π ν n tanδ / c          (Lau 2006: tanδ power loss per radian)
//   ABSLENGTH(ν) = 1/α = c / (2π ν n tanδ)
// so ABSLENGTH ∝ 1/ν. Tabulated on the same 10 GHz–20 THz log grid as the
// Cu REFLECTIVITY / Planck CDF.
// ---------------------------------------------------------------------------
inline G4Material* BuildDielectricMaterial(const G4String& name,
                                            const G4String& nistBase,
                                            G4double n_index,
                                            G4double tan_delta)
{
  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat) return mat;

  // Clone the NIST base under our own name so distinct (n, tanδ) tables never
  // collide on the shared NIST instance, and the ROOT legend records a
  // meaningful material name.
  mat = G4NistManager::Instance()->BuildMaterialWithNewDensity(name, nistBase);

  const G4double c_SI  = 2.998e8;             // m/s
  const G4double h_eVs = 4.13566769692e-15;   // eV·s

  const int      N    = 24;
  const G4double Emin = 4.14e-5*eV;   // 10 GHz
  const G4double Emax = 8.27e-2*eV;   // 20 THz

  std::vector<G4double> energies(N), abslen(N);
  const G4double logMin = std::log(Emin), logMax = std::log(Emax);
  for (int i = 0; i < N; ++i) {
    energies[i]         = std::exp(logMin + i*(logMax - logMin)/(N - 1.));
    const G4double nu   = (energies[i]/eV) / h_eVs;                 // Hz
    const G4double alpha = 2.*CLHEP::pi * nu * n_index * tan_delta / c_SI; // 1/m
    abslen[i]           = (1./alpha) * m;     // 1/m → Geant4 length units
  }

  const std::vector<G4double> e2 = {Emin, Emax};
  const std::vector<G4double> ri = {n_index, n_index};

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",    e2,       ri);
  mpt->AddProperty("ABSLENGTH", energies, abslen);
  mat->SetMaterialPropertiesTable(mpt);
  return mat;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

// Near-vacuum material that flags crack volumes for HFSS diffraction routing.
// RINDEX=1, no REFLECTIVITY — identified by name "vacuum_wg" in PostStepDoIt.
inline G4Material* GetVacuumWG()
{
  G4Material* mat = G4Material::GetMaterial("vacuum_wg", false);
  if (mat) return mat;
  mat = new G4Material("vacuum_wg", 1., 1.008*g/mole, 1e-25*g/cm3);
  const std::vector<G4double> e  = {1e-6*eV, 1.0*eV};
  const std::vector<G4double> ri = {1.0, 1.0};
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

// Copper with full Drude reflectance model, parameterized by (RRR, T_K).
// Material name is deterministic ("Cu_RRR{RRR}_T{T_K}K", %g formatting so
// e.g. 4 K → "T4K", 4.6 K → "T4.6K") so repeated calls with the same
// arguments return the cached G4Material and distinct temperatures never
// alias to the same cache entry.
// T_K defaults to 4.0 K (cryogenic baseline for BBRsim).
// Below 50 K: σ_DC = RRR × σ_RT (phonons frozen).  Above 50 K: Matthiessen.
inline G4Material* GetCopper(G4int RRR, G4double T_K = 4.0)
{
  char tbuf[32];
  std::snprintf(tbuf, sizeof(tbuf), "%g", T_K);
  G4String name = "Cu_RRR" + std::to_string(RRR) + "_T" + tbuf + "K";
  return BuildDrudeMaterial(name, 29., 63.546, 8.96, RRR, T_K);
}

// Name-to-RRR mapping for backward compatibility and messenger use.
// RRR values derived from Serov 2016 measured D at 4 K:
//   OFHC_Cu: σ_eff = 5.96×10⁹ S/m → RRR = 100
//   OF_Cu:   σ_eff = 1.79×10⁸ S/m → RRR =   3
//   HP_Cu:   σ_eff = 3.38×10⁸ S/m → RRR =   6
inline G4Material* GetCopperByName(const G4String& name)
{
  if (name == "OFHC_Cu") return GetCopper(100, 4.0);
  if (name == "OF_Cu")   return GetCopper(  3, 4.0);
  if (name == "HP_Cu")   return GetCopper(  6, 4.0);
  return nullptr;
}

// Cirlex (pressure-formed DuPont Kapton polyimide) at cryogenic temperature.
// n ≈ 1.95, tan δ ≈ 0.015 — single complex-permittivity fit at 5 K over
// 300 GHz–3 THz (Lau et al. 2006, Appl. Opt. 45, 3746; corroborated by the
// CMBPol optics review, Table 1). Base material G4_KAPTON.
inline G4Material* GetCirlex()
{
  return BuildDielectricMaterial("Cirlex", "G4_KAPTON", 1.95, 0.015);
}

// Crystalline silicon detector substrate. n = 3.39 (flat in the trans-mm
// band, YYC / Frey NASA Goddard); tan δ = 1e-4 (Chang §5.3.1.4).
inline G4Material* GetSiliconCrystal()
{
  return BuildDielectricMaterial("Si", "G4_Si", 3.39, 1.0e-4);
}

// Crystalline germanium detector substrate. n ≈ 4.0 (flat in the trans-mm
// band); tan δ = 6e-5 (Chang §5.3.1.4).
inline G4Material* GetGermaniumCrystal()
{
  return BuildDielectricMaterial("Ge", "G4_Ge", 4.0, 6.0e-5);
}

} // namespace BBRMaterials

#endif
