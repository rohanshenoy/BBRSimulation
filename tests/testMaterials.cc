// testMaterials.cc — BBRMaterials: the Cu REFLECTIVITY table against the
// complex Drude formula (re-implemented here with the header's constants, so
// it catches drift) and against physics anchors that do not share the code
// (the relaxation plateau, Hagen-Rubens, Serov 2016); the dielectric
// ABSLENGTH tables, on and between the grid points; the flag materials; the
// input guards BBR015 / BBR016.
#include "BBRTestSupport.hh"

#include "BBRMaterials.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <complex>
#include <iterator>
#include <vector>

namespace {

// The header's constants (BBRMaterials.hh).
const double kSigmaRT = 5.96e7;              // S/m
const double kNe      = 8.49e28;             // m^-3
const double kMe      = 9.109e-31;           // kg
const double kE       = 1.602e-19;           // C
const double kEps0    = 8.8541878128e-12;    // F/m
const double kH_eVs   = 4.13566769692e-15;   // eV s
const double kPi      = 3.14159265358979323846;
const int    kRRRs[]  = {1, 3, 6, 100};

G4MaterialPropertyVector* Refl(G4Material* m) {
  return m->GetMaterialPropertiesTable()->GetProperty("REFLECTIVITY");
}

double Tau(int RRR) { return RRR * kSigmaRT * kMe / (kNe * kE * kE); }

// Absorptance D = 1 - R of the full complex Drude model below 50 K
// (sigma_DC = RRR sigma_RT), Griffiths §9.4 with complex sigma. Written as
// D = 4 Re(n) / |n + 1|^2, which equals 1 - |(n-1)/(n+1)|^2 without the
// cancellation of 1 - R near R = 1.
double DrudeD(int RRR, double nu_Hz) {
  const double w = 2. * kPi * nu_Hz;
  const std::complex<double> sigma =
    RRR * kSigmaRT / std::complex<double>(1., -w * Tau(RRR));
  const std::complex<double> eps = 1. + std::complex<double>(0., 1.) * sigma / (kEps0 * w);
  const std::complex<double> n = std::sqrt(eps);
  return 4. * n.real() / std::norm(n + 1.);
}

double Nu(double E) { return (E / eV) / kH_eVs; }
double E_of(double nu_Hz) { return kH_eVs * nu_Hz * eV; }

struct Dielectric {
  const char* name;
  G4Material* (*get)();
  double n, tand;
  const char* base;
  double grid0_mm;
};

// The three loss-tangent dielectrics as BBRMaterials.hh builds them, with the
// NIST base and the ABSLENGTH at grid point 0.
const Dielectric kDielectrics[] = {
  {"Si", BBRMaterials::GetSiliconCrystal, 3.39, 1.0e-4, "G4_Si", 1.4060389658e4},
  {"Ge", BBRMaterials::GetGermaniumCrystal, 4.0, 6.0e-5, "G4_Ge", 1.9860300392e4},
  {"Cirlex", BBRMaterials::GetCirlex, 1.95, 0.015, "G4_KAPTON", 162.95631091},
};

// ABSLENGTH = c / (2 pi nu n tan delta) in mm at photon energy E, with the
// header's c = 2.998e8 m/s.
double AbsLen_mm(const Dielectric& d, double E) {
  return 2.998e8 / (2. * kPi * Nu(E) * d.n * d.tand) * 1e3;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace bbrtest;
  return RunCase(argc, argv, {
    {"cu_table_shape", [] {
      for (int rrr : kRRRs) {
        auto* R = Refl(BBRMaterials::GetCopper(rrr, 4.));
        CHECK(R->GetVectorLength() == 24);
        CHECK_REL(R->Energy(0) / eV, 4.14e-5, 1e-12);
        CHECK_REL(R->Energy(23) / eV, 8.27e-2, 1e-12);
        double maxdev = 0.;
        for (std::size_t i = 0; i + 1 < R->GetVectorLength(); ++i)
          maxdev = std::max(maxdev, std::abs(R->Energy(i + 1) / R->Energy(i) - 1.391554514868290));
        CHECK_NEAR(maxdev, 0., 1e-12);
      }
    }},
    {"cu_drude_grid", [] {
      const double d0[]  = {2.731580e-4, 1.574700e-4, 1.110892e-4, 2.528803e-5};
      const double d23[] = {4.812452e-3, 1.624294e-3, 8.133363e-4, 4.883606e-5};
      int k = 0;
      for (int rrr : kRRRs) {
        auto* R = Refl(BBRMaterials::GetCopper(rrr, 4.));
        double worst = 0.;
        for (std::size_t i = 0; i < R->GetVectorLength(); ++i) {
          const double D = 1. - (*R)[i], Dx = DrudeD(rrr, Nu(R->Energy(i)));
          worst = std::max(worst, std::abs(D - Dx) / Dx);
        }
        CHECK_NEAR(worst, 0., 1e-9);
        CHECK_REL(1. - (*R)[0], d0[k], 1e-5);
        CHECK_REL(1. - (*R)[23], d23[k], 1e-5);
        ++k;
      }
    }},
    {"cu_limits", [] {
      // Relaxation plateau D = 2 / (omega_p tau) for RRR 100 at 20 THz.
      const double wp = std::sqrt(kNe * kE * kE / (kEps0 * kMe));
      CHECK_REL(wp, 1.6436e16, 1e-4);
      auto* R100 = Refl(BBRMaterials::GetCopper(100, 4.));
      CHECK_REL(1. - (*R100)[23], 2. / (wp * Tau(100)), 1e-3);
      // Hagen-Rubens D = 2 sqrt(2 eps0 omega / sigma) for RRR 1 at 10 GHz.
      auto* R1 = Refl(BBRMaterials::GetCopper(1, 4.));
      const double w0 = 2. * kPi * Nu(R1->Energy(0));
      CHECK_REL(1. - (*R1)[0], 2. * std::sqrt(2. * kEps0 * w0 / kSigmaRT), 2e-3);
      // D strictly decreasing with RRR at every grid point.
      std::vector<G4MaterialPropertyVector*> Rs;
      for (int rrr : kRRRs) Rs.push_back(Refl(BBRMaterials::GetCopper(rrr, 4.)));
      bool decreasing = true;
      for (std::size_t i = 0; i < 24; ++i)
        for (std::size_t j = 0; j + 1 < Rs.size(); ++j)
          decreasing = decreasing && (1. - (*Rs[j])[i]) > (1. - (*Rs[j + 1])[i]);
      CHECK(decreasing);
    }},
    {"cu_serov_of", [] {
      // Serov 2016: OF Cu (RRR 3) at 150 GHz, D = 0.58e-3. HP_Cu (RRR 6 vs
      // Serov at 230 GHz) is the runner's XFAIL and is not pinned here.
      const double D = 1. - Refl(BBRMaterials::GetCopper(3, 4.))->Value(E_of(150e9));
      CHECK_REL(D, 5.80e-4, 0.10);
    }},
    {"cu_offgrid_interp", [] {
      // Geant4 interpolates linearly in energy, so off-grid D is up to 0.52 %
      // off the formula (a log-log interpolation would give 0.26 %).
      auto* R100 = Refl(BBRMaterials::GetCopper(100, 4.));
      CHECK_REL(DrudeD(100, 500e9), 4.873578e-5, 1e-6);
      CHECK_REL(1. - R100->Value(E_of(500e9)), 4.873578e-5, 6e-3);
      for (int rrr : kRRRs) {
        auto* R = Refl(BBRMaterials::GetCopper(rrr, 4.));
        const double l0 = std::log(R->Energy(0)), l1 = std::log(R->Energy(23));
        double worst = 0.;
        for (int i = 0; i < 1000; ++i) {
          const double E = std::exp(l0 + (l1 - l0) * (i + 0.5) / 1000.);
          const double Dx = DrudeD(rrr, Nu(E));
          worst = std::max(worst, std::abs((1. - R->Value(E)) - Dx) / Dx);
        }
        CHECK(worst <= 5.2e-3);
        std::printf("RRR %d: worst off-grid |dD/D| = %.4e\n", rrr, worst);
        CHECK(R->Value(R->Energy(0) / 2.) == (*R)[0]);
        CHECK(R->Value(R->Energy(23) * 2.) == (*R)[23]);
      }
    }},
    {"cu_cache", [] {
      G4Material* a = BBRMaterials::GetCopper(100, 4.);
      CHECK(BBRMaterials::GetCopper(100, 4.) == a);
      G4Material* b = BBRMaterials::GetCopper(100, 4.6);
      CHECK(b != a);
      CHECK(b->GetName() == "Cu_RRR100_T4.6K");
    }},
    {"cu_invalid_input", [] {  // D8
      ExpectG4Exception("BBR015", [] { BBRMaterials::GetCopper(0, 4.); }, "BuildDrudeMaterial");
      ExpectG4Exception("BBR015", [] { BBRMaterials::GetCopper(-5, 4.); }, "BuildDrudeMaterial");
      ExpectG4Exception("BBR015", [] { BBRMaterials::GetCopper(100, 0.); }, "BuildDrudeMaterial");
      ExpectG4Exception("BBR015", [] { BBRMaterials::GetCopper(100, std::nan("")); }, "BuildDrudeMaterial");
    }},
    {"dielectrics", [] {
      for (const auto& d : kDielectrics) {
        G4Material* m = d.get();
        CHECK(m->GetName() == d.name);
        auto* mpt = m->GetMaterialPropertiesTable();
        auto* A = mpt->GetProperty("ABSLENGTH");
        auto* n = mpt->GetProperty("RINDEX");
        CHECK(A && A->GetVectorLength() == 24);
        double worst = 0.;
        for (std::size_t i = 0; i < A->GetVectorLength(); ++i) {
          const double ref = AbsLen_mm(d, A->Energy(i));
          worst = std::max(worst, std::abs((*A)[i] / mm - ref) / ref);
        }
        CHECK_NEAR(worst, 0., 1e-12);
        CHECK_REL((*A)[0] / mm, d.grid0_mm, 1e-8);
        CHECK(n && n->GetVectorLength() == 2 && (*n)[0] == d.n && (*n)[1] == d.n);
        CHECK(mpt->GetProperty("REFLECTIVITY") == nullptr);
        CHECK_REL(m->GetDensity(),
                  G4NistManager::Instance()->FindOrBuildMaterial(d.base)->GetDensity(), 1e-12);
      }
    }},
    {"dielectric_offgrid_interp", [] {
      // Geant4 interpolates ABSLENGTH linearly in energy, and c/(2 pi nu n tand)
      // is convex in E, so between grid points E_i and r E_i the table lies above
      // the formula, by at most (r + 1)^2 / (4 r) - 1 = 2.754 % at the arithmetic
      // midpoint (2.736 % at the log midpoint). RINDEX is flat, so it is exact.
      const double r = 1.391554514868290, bound = (r + 1.) * (r + 1.) / (4. * r) - 1.;
      CHECK_REL(bound, 2.75438e-2, 1e-5);
      // The formula and the table at 500 GHz [mm], in kDielectrics order.
      struct At500 { double formula_mm, table_mm; };
      const At500 at500[] = {{281.50237, 286.21944}, {397.62210, 404.28497}, {3.2625403, 3.3172100}};
      static_assert(std::size(at500) == std::size(kDielectrics));
      for (std::size_t k = 0; k < std::size(kDielectrics); ++k) {
        const Dielectric& c = kDielectrics[k];
        auto* mpt = c.get()->GetMaterialPropertiesTable();
        auto* A = mpt->GetProperty("ABSLENGTH");
        auto* n = mpt->GetProperty("RINDEX");
        auto L = [&](double E) { return AbsLen_mm(c, E); };
        // Arithmetic midpoint of every interval: exactly (r_i + 1)^2 / (4 r_i).
        double worstMid = 0.;
        for (std::size_t i = 0; i + 1 < A->GetVectorLength(); ++i) {
          const double E0 = A->Energy(i), E1 = A->Energy(i + 1), ri = E1 / E0, Em = 0.5 * (E0 + E1);
          worstMid = std::max(worstMid, std::abs(A->Value(Em) / mm / L(Em) / ((ri + 1.) * (ri + 1.) / (4. * ri)) - 1.));
        }
        CHECK_NEAR(worstMid, 0., 1e-12);
        // 1000 log-spaced energies inside the grid: never below the formula,
        // never more than the bound above it, and within 1e-4 of the bound at worst.
        const double l0 = std::log(A->Energy(0)), l1 = std::log(A->Energy(23));
        double lo = 1., hi = -1.;
        bool flatN = true;
        for (int i = 0; i < 1000; ++i) {
          const double E = std::exp(l0 + (l1 - l0) * (i + 0.5) / 1000.);
          const double x = A->Value(E) / mm / L(E) - 1.;
          lo = std::min(lo, x);
          hi = std::max(hi, x);
          flatN = flatN && n->Value(E) == c.n;
        }
        std::printf("%s: off-grid ABSLENGTH / formula - 1 in [%.3e, %.5e]\n", c.name, lo, hi);
        CHECK(lo >= -1e-12);
        CHECK(hi <= bound + 1e-12);
        CHECK(hi >= bound - 1e-4);
        CHECK(flatN);
        // 500 GHz lies between grid points 11 and 12: +1.676 % for each of them.
        const double E500 = E_of(500e9);
        CHECK_REL(L(E500), at500[k].formula_mm, 1e-7);
        CHECK_REL(A->Value(E500) / mm, at500[k].table_mm, 1e-7);
        CHECK_REL(A->Value(E500) / mm / L(E500) - 1., 1.67568e-2, 1e-4);
        // Outside the grid: the edge values.
        CHECK(A->Value(A->Energy(0) / 2.) == (*A)[0]);
        CHECK(A->Value(A->Energy(23) * 2.) == (*A)[23]);
        CHECK(n->Value(A->Energy(0) / 2.) == c.n && n->Value(A->Energy(23) * 2.) == c.n);
      }
    }},
    {"dielectric_name_collision", [] {  // D7: a user-made "Si" without optics
      new G4Material("Si", 14., 28.0855 * g / mole, 2.33 * g / cm3);
      ExpectG4Exception("BBR016", [] { BBRMaterials::GetSiliconCrystal(); },
                        "BuildDielectricMaterial");
    }},
    {"flag_materials", [] {
      struct Flag { G4Material* (*get)(); const char* name; int refl; };  // refl: -1 none
      const Flag fs[] = {{BBRMaterials::GetVacuumWG, "vacuum_wg", -1},
                         {BBRMaterials::GetPerfectAbsorber, "BBR_PerfectAbsorber", 0},
                         {BBRMaterials::GetPerfectReflector, "BBR_PerfectReflector", 1}};
      for (const auto& f : fs) {
        G4Material* m = f.get();
        CHECK(m->GetName() == f.name);
        CHECK(f.get() == m);
        auto* mpt = m->GetMaterialPropertiesTable();
        auto* ri = mpt->GetProperty("RINDEX");
        CHECK(ri && ri->GetVectorLength() == 2);
        CHECK_REL(ri->Energy(0) / eV, 1e-6, 1e-12);
        CHECK_REL(ri->Energy(1) / eV, 1., 1e-12);
        CHECK((*ri)[0] == 1. && (*ri)[1] == 1.);
        auto* r = mpt->GetProperty("REFLECTIVITY");
        if (f.refl < 0) {
          CHECK(r == nullptr);
        } else {
          CHECK(r != nullptr);
          if (r) CHECK(r->Value(E_of(20e12)) == double(f.refl));
        }
      }
    }},
  });
}
