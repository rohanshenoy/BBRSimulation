// testWrapper — BBSimOpBoundaryProcess::PostStepDoIt on hand-built, navigated
// steps (tests/WrapperWorld.hh): the HFSS diffraction path. Expected values
// are derived by hand from the fixture tables (tests/HFSSFixture.hh) and the
// HFSS frame convention: crack axis +x, long axis +y (theta_hat), gap +z
// (phi_hat); IWaveTheta = acos(-k.n), IWavePhi = atan2(-k.theta_hat, k.phi_hat)
// folded into [0, 90] by the mirror signs sy (theta_hat) and sx (phi_hat).
#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "WrapperWorld.hh"

#include "BBRConfigManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4coutDestination.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "CLHEP/Random/MixMaxRng.h"

#include <cfloat>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace bbrtest;
using wrapperworld::EnergyEV;
using wrapperworld::World;
using S = BBSimOpBoundaryProcess;

namespace {
const G4ThreeVector X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1);
const double c45 = std::cos(45 * deg), s45 = std::sin(45 * deg);
const G4ThreeVector kHitMinusX(-2 * mm, 0, 0);   // centre of the crack's -x face
const double E500 = EnergyEV(500.);

double Inset() { return 20. * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(); }

// Writes the dataset as <root>/waveguides/crack_<freq>GHz and points the library at root.
void UseData(TempDir& d, const hfssfix::Dataset& ds, const std::string& freq = "500") {
  hfssfix::WriteDataset(d.path(), "crack_" + freq + "GHz", ds);
  BBRConfigManager::SetDataDir(d.path().string());
}

// Replaces whatever is left in the engine by the given draws.
void Script(ScriptedEngine& e, std::initializer_list<double> u) {
  while (e.Remaining()) e.flat();
  for (double v : u) e.Push(v);
}

// Polarization in the plane of incidence (the x-k plane), perpendicular to k.
G4ThreeVector InPlane(const G4ThreeVector& k) { return (X - X.dot(k.unit()) * k.unit()).unit(); }

// A slab material carrying a material REFLECTIVITY table (what BBRMaterials::GetCopper builds).
G4Material* Mirror(const char* name, const std::vector<G4double>& e, const std::vector<G4double>& r) {
  auto* m = new G4Material(name, 29., 63.546 * g / mole, 8.96 * g / cm3);
  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("REFLECTIVITY", e, r);
  m->SetMaterialPropertiesTable(mpt);
  return m;
}
G4Material* Flat30() { return Mirror("R30", {1e-6 * eV, 1 * eV}, {0.3, 0.3}); }
const G4ThreeVector k30(std::cos(30 * deg), std::sin(30 * deg), 0);
const G4ThreeVector kSlabStart(-10 * mm, -15 * mm, 0);   // k30 from here hits the slab face x = -2 mm at y = -10.4 mm

// Direction of far-field row (Phi, Theta) in the folded frame: sinT cosP n + sinT sinP tf + cosT pf.
G4ThreeVector FFDir(double P, double T, const G4ThreeVector& n, const G4ThreeVector& tf, const G4ThreeVector& pf) {
  P *= deg; T *= deg;
  return std::sin(T) * std::cos(P) * n + std::sin(T) * std::sin(P) * tf + std::cos(T) * pf;
}

// --- exit start (CDX-027) -----------------------------------------------------
// The start on the ray through the sampled exit point pos_out (on the exit plane, whose outward
// axis is nOut) at the normal inset: its normal distance from the plane is exactly Inset().
G4ThreeVector ExpectedStart(const G4ThreeVector& pos_out, const G4ThreeVector& dir, const G4ThreeVector& nOut) {
  return pos_out - (Inset() / dir.dot(nOut)) * dir;
}

// True if the straight segment from r.pos along r.dir first leaves the box of half-extents half
// (centred at the origin, unrotated) through the plane x = exitX, with |y| < half.y and |z| < half.z
// there; the crossing point goes to *crossing.
bool LeavesThroughExitFace(const wrapperworld::Result& r, const G4ThreeVector& half, double exitX,
                           G4ThreeVector* crossing = nullptr) {
  const G4ThreeVector& p = r.pos;
  const G4ThreeVector& v = r.dir;
  if (std::abs(p.x()) >= half.x() || std::abs(p.y()) >= half.y() || std::abs(p.z()) >= half.z()) return false;
  auto tFace = [](double pc, double vc, double h) { return vc > 0 ? (h - pc) / vc : vc < 0 ? (-h - pc) / vc : DBL_MAX; };
  const double tx = tFace(p.x(), v.x(), half.x()), ty = tFace(p.y(), v.y(), half.y()), tz = tFace(p.z(), v.z(), half.z());
  if (!(tx < ty && tx < tz)) return false;
  const G4ThreeVector c = p + tx * v;
  if (crossing) *crossing = c;
  return std::abs(c.x() - exitX) < 1e-12 * mm && std::abs(c.y()) < half.y() && std::abs(c.z()) < half.z();
}

// The same for a cylinder of radius R about the x axis, half-length halfX.
bool LeavesThroughDiscExitFace(const wrapperworld::Result& r, double R, double halfX, double exitX,
                               G4ThreeVector* crossing = nullptr) {
  const G4ThreeVector& p = r.pos;
  const G4ThreeVector& v = r.dir;
  if (std::abs(p.x()) >= halfX || std::hypot(p.y(), p.z()) >= R) return false;
  const double tx = v.x() > 0 ? (halfX - p.x()) / v.x() : v.x() < 0 ? (-halfX - p.x()) / v.x() : DBL_MAX;
  const double a = v.y() * v.y() + v.z() * v.z(), b = 2 * (p.y() * v.y() + p.z() * v.z()),
               c = p.y() * p.y() + p.z() * p.z() - R * R;
  const double tr = a > 0 ? (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a) : DBL_MAX;
  if (!(tx < tr)) return false;
  const G4ThreeVector x = p + tx * v;
  if (crossing) *crossing = x;
  return std::abs(x.x() - exitX) < 1e-12 * mm && std::hypot(x.y(), x.z()) < R;
}

// Writes one dataset per frequency label as crack_<label>GHz (relabelled from 500 GHz) and points the
// library at root. discR > 0: the sidecars declare a disc of that radius (m), outside points omitted.
void UseGrid(TempDir& d, const std::vector<std::pair<std::string, hfssfix::Dataset>>& sets, double discR = 0.) {
  for (const auto& [label, ds0] : sets) {
    const std::string stem = "crack_" + label + "GHz";
    const hfssfix::Dataset ds = hfssfix::Relabel(ds0, label);
    const auto wg = hfssfix::WriteDataset(d.path(), stem, ds, discR <= 0.);
    if (discR > 0.) {
      auto p = hfssfix::SidecarFrom(stem, ds);
      p.disc = true; p.radius = discR; p.outsidePoints = "omitted";
      WriteFile(wg / (stem + ".dataset.json"), hfssfix::SidecarJson(p));
    }
  }
  BBRConfigManager::SetDataDir(d.path().string());
}

// Far-field (Phi, Theta) [deg] of the outgoing directions n X + s Z (toward the +z wall), n X - s Z
// (away from it) and n X + s Y (toward the +y wall), s = sqrt(1 - n^2).
double Deg(double rad) { return rad * 180. / CLHEP::pi; }
std::pair<double, double> TowardG(double n) { return {0., Deg(std::asin(n))}; }
std::pair<double, double> AwayG(double n) { return {0., 180. - Deg(std::asin(n))}; }
std::pair<double, double> TowardL(double n) { return {Deg(std::acos(n)), 90.}; }
hfssfix::Dataset OnePoint(double y, double z, std::pair<double, double> pt) {
  return hfssfix::OnePointOneDirection(y, z, pt.first, pt.second);
}

// The crack equal to the fixtures' declared section (OR-020), so a clamped wall row lies 1 nm from its wall.
const G4ThreeVector kExact(2 * mm, 4.5 * mm, 0.025 * mm);
const double kYw = (4.5e-3 - 1e-9) * CLHEP::m, kZw = (2.5e-5 - 1e-9) * CLHEP::m;   // clamped wall rows

// Collects what G4cout receives while installed.
struct Capture : public G4coutDestination {
  std::string out;
  G4int ReceiveG4cout(const G4String& s) override { out += s; return 0; }
};
}  // namespace

int main(int argc, char** argv) {
  return RunCase(argc, argv, {
    {"fold_six_azimuths", [] {
      TempDir d; UseData(d, hfssfix::WrapKeys());
      World w; ScriptedEngine eng;
      struct Case { G4ThreeVector k; double sy, sx, P, T; };
      // kt = -k_y, kp = k_z; phi_raw = atan2(kt, kp); far-field row of the selected key
      const Case cases[] = {
        {{c45, -s45, 0}, +1, +1, 60, 90},     // phi_raw = +90  -> (90,135)
        {{c45, +s45, 0}, -1, +1, 60, 90},     // phi_raw = -90  -> (90,135)
        {{c45, 0, +s45}, +1, +1, 30, 60},     // phi_raw = 0    -> (0,135)
        {{c45, 0, -s45}, +1, -1, 30, 60},     // phi_raw = 180  -> (0,135)
        {{c45, -0.5, 0.5}, +1, +1, 0, 60},    // phi_raw = +45  -> (45,135)
        {{c45, 0.5, -0.5}, -1, -1, 0, 60},    // phi_raw = -135 -> (45,135)
      };
      for (const auto& c : cases) {
        Script(eng, {0.0, 0.5, 0.5});   // transmit, direction, position
        const auto r = w.Shoot(kHitMinusX - 1 * mm * c.k, c.k, InPlane(c.k), E500, eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(2 * mm, c.sy * 1 * mm, c.sx * 0.01 * mm),
                                       FFDir(c.P, c.T, X, c.sy * Y, c.sx * Z), X), 1e-9);
        CHECK_VEC(r.dir, FFDir(c.P, c.T, X, c.sy * Y, c.sx * Z), 1e-9);
      }
    }},
    {"exit_start_wall_grazing", [] {
      // A clamped +z wall row (z = 25 um - 1 nm) leaving at a grazing angle toward that wall: from the
      // normal inset the segment would meet z = +25 um after ~1 nm, before the exit plane (CDX-027). The
      // start on the ray through the sampled point leaves through the exit face exactly at that point.
      TempDir d;
      const double ns[] = {0.02, 1e-2, 1e-3};
      const char* labels[] = {"100", "200", "400"};
      const double freqs[] = {100., 200., 400.};
      UseGrid(d, {{labels[0], OnePoint(0., 2.5e-5, TowardG(ns[0]))},
                  {labels[1], OnePoint(0., 2.5e-5, TowardG(ns[1]))},
                  {labels[2], OnePoint(0., 2.5e-5, TowardG(ns[2]))}});
      World w(kExact); ScriptedEngine eng;
      for (int i = 0; i < 3; ++i) {
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, EnergyEV(freqs[i]), eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_NEAR(r.dir.x(), ns[i], 1e-12);
        CHECK(r.dir.z() > 0.);
        G4ThreeVector c;
        CHECK(LeavesThroughExitFace(r, kExact, 2 * mm, &c));
        CHECK_VEC(c, G4ThreeVector(2 * mm, 0, kZw), 1e-12 * mm);
        CHECK_VEC(r.pos, c - (Inset() / r.dir.x()) * r.dir, 1e-12 * mm);
      }
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
    }},
    {"exit_start_away_from_wall", [] {
      // The same wall row leaving grazing AWAY from the wall: the back-off along the ray would cross the
      // wall behind it, so the start is the normal inset, whose ray heads into the crack.
      TempDir d; UseGrid(d, {{"500", OnePoint(0., 2.5e-5, AwayG(1e-3))}});
      World w(kExact); ScriptedEngine eng;
      Script(eng, {0.0, 0.5, 0.5});
      const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit);
      CHECK(r.dir.z() < 0.);
      CHECK_VEC(r.pos, G4ThreeVector(2 * mm - Inset(), 0, kZw), 1e-12 * mm);
      CHECK(LeavesThroughExitFace(r, kExact, 2 * mm));
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
      CHECK(w.Process().GetExitStartNormalCount() == 1);
    }},
    {"exit_start_long_wall", [] {
      // A clamped +y wall row (y = 4.5 mm - 1 nm) leaving grazing toward that wall (T = 90, P = acos n).
      TempDir d; UseGrid(d, {{"500", OnePoint(4.5e-3, 0., TowardL(1e-3))}});
      World w(kExact); ScriptedEngine eng;
      Script(eng, {0.0, 0.5, 0.5});
      const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit);
      CHECK(r.dir.y() > 0.);
      G4ThreeVector c;
      CHECK(LeavesThroughExitFace(r, kExact, 2 * mm, &c));
      CHECK_VEC(c, G4ThreeVector(2 * mm, kYw, 0), 1e-12 * mm);
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
    }},
    {"exit_start_folds", [] {
      // fold_six_azimuths' incidences: the mirrors (sy, sx) carry the wall row and the wall-ward direction
      // to the same wall, which the photon must still leave beside, through the exit face.
      TempDir d; UseGrid(d, {{"500", OnePoint(0., 2.5e-5, TowardG(1e-3))}});
      World w(kExact); ScriptedEngine eng;
      struct Case { G4ThreeVector k; double sy, sx; };
      const Case cases[] = {
        {{c45, -s45, 0}, +1, +1}, {{c45, +s45, 0}, -1, +1}, {{c45, 0, +s45}, +1, +1},
        {{c45, 0, -s45}, +1, -1}, {{c45, -0.5, 0.5}, +1, +1}, {{c45, 0.5, -0.5}, -1, -1},
      };
      for (const auto& c : cases) {
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot(kHitMinusX - 1 * mm * c.k, c.k, InPlane(c.k), E500, eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.dir, FFDir(TowardG(1e-3).first, TowardG(1e-3).second, X, c.sy * Y, c.sx * Z), 1e-12);
        CHECK(LeavesThroughExitFace(r, kExact, 2 * mm));
        CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(2 * mm, 0, c.sx * kZw), r.dir, X), 1e-12 * mm);
      }
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
    }},
    {"exit_start_minus_x_face", [] {
      // Entry through the +x face: the exit face is x = -2 mm and the same rule holds, mirrored.
      TempDir d; UseGrid(d, {{"500", OnePoint(0., 2.5e-5, TowardG(1e-3))}});
      World w(kExact); ScriptedEngine eng;
      Script(eng, {0.0, 0.5, 0.5});
      const auto r = w.Shoot({3 * mm, 0, 0}, -X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit);
      CHECK_NEAR(r.dir.x(), -1e-3, 1e-12);
      G4ThreeVector c;
      CHECK(LeavesThroughExitFace(r, kExact, -2 * mm, &c));
      CHECK_VEC(c, G4ThreeVector(-2 * mm, 0, kZw), 1e-12 * mm);
      CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(-2 * mm, 0, kZw), r.dir, -X), 1e-12 * mm);
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
    }},
    {"exit_start_guard_limit", [] {
      // n = 2e-6, just above the samplers' 1e-6 outward guard, toward the +z wall: no exception, and the
      // photon either leaves through the exit face or is counted as one fallback.
      TempDir d; UseGrid(d, {{"500", OnePoint(0., 2.5e-5, TowardG(2e-6))}});
      World w(kExact); ScriptedEngine eng;
      Script(eng, {0.0, 0.5, 0.5});
      const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit);
      CHECK_NEAR(r.dir.x(), 2e-6, 1e-12);
      const G4int fb = w.Process().GetExitStartFallbackCount();
      CHECK(fb == 0 || fb == 1);
      CHECK(LeavesThroughExitFace(r, kExact, 2 * mm) || fb == 1);
    }},
    {"exit_start_disc_rim", [] {
      // The round gap (R = 50 um): the rim node (0, R), clamped to R - 1 nm, leaving at n = 1e-3 radially
      // outward (B valid; from the normal inset the segment meets the rim) and tangentially (B valid).
      const double R = 5e-5;
      TempDir d;
      UseGrid(d, {{"100", hfssfix::OnePointOneDirectionDisc(R, 0., R, TowardG(1e-3).first, TowardG(1e-3).second)},
                  {"200", hfssfix::OnePointOneDirectionDisc(R, 0., R, TowardL(1e-3).first, TowardL(1e-3).second)}},
             R);
      World w({2 * mm, 1 * mm, 1 * mm}, {}, nullptr, wrapperworld::RoundGapSolid(R * CLHEP::m, 2 * mm));
      ScriptedEngine eng;
      const G4ThreeVector rim(2 * mm, 0, (R - 1e-9) * CLHEP::m);
      for (double f : {100., 200.}) {
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, EnergyEV(f), eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_NEAR(r.dir.x(), 1e-3, 1e-12);
        CHECK(f == 100. ? r.dir.z() > 0.9 : r.dir.y() > 0.9);
        G4ThreeVector c;
        CHECK(LeavesThroughDiscExitFace(r, R * CLHEP::m, 2 * mm, 2 * mm, &c));
        CHECK_VEC(c, rim, 1e-12 * mm);
        CHECK_VEC(r.pos, ExpectedStart(rim, r.dir, X), 1e-12 * mm);
      }
      CHECK(w.Process().GetExitStartFallbackCount() == 0);
    }},
    {"exit_start_corner", [] {
      // The corner node (+y wall, +z wall) leaving toward +z and away from +y at n = 1e-3: the ray start
      // crosses the +y wall and the normal-inset ray meets the +z wall, so the photon keeps the normal inset
      // and is counted; the count is printed on the diffraction summary line, every hundredth event.
      const double n = 1e-3, a = std::sqrt(1. - n * n) / std::sqrt(2.);
      const double T = Deg(std::acos(a)), P = Deg(std::atan2(-a, n));
      TempDir d; UseGrid(d, {{"500", hfssfix::OnePointOneDirection(4.5e-3, 2.5e-5, P, T)}});
      World w(kExact); ScriptedEngine eng;
      Script(eng, {0.0, 0.5, 0.5});
      const auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit);
      CHECK_VEC(r.dir, G4ThreeVector(n, -a, a), 1e-12);
      CHECK(w.Process().GetExitStartFallbackCount() == 1);
      CHECK(w.Process().GetExitStartNormalCount() == 0);
      CHECK_VEC(r.pos, G4ThreeVector(2 * mm - Inset(), kYw, kZw), 1e-12 * mm);
      Capture cap;
      G4iosSetDestination(&cap);
      for (int i = 0; i < 99; ++i) {
        Script(eng, {0.0, 0.5, 0.5});
        w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      }
      G4iosSetDestination(nullptr);
      CHECK(w.Process().GetExitStartFallbackCount() == 100);
      CHECK(w.Process().GetExitStartNormalCount() == 0);
      CHECK(cap.out.find("[BBR] diffraction events=100 T_obs=1 exit_start_fallback=99 exit_start_normal=0\n") !=
            std::string::npos);
    }},
    {"run_totals", [] {
      // Process-wide totals over all threads: HFSS transmissions, exit-start fallbacks and starts at the
      // normal inset with a valid ray, reset by the run. One grid holds the three cases, one frequency each
      // (the crack library caches a dataset by ID for the life of the process).
      using Totals = BBSimOpBoundaryProcess::Totals;
      const double n = 1e-3, a = std::sqrt(1. - n * n) / std::sqrt(2.);
      const double T = Deg(std::acos(a)), P = Deg(std::atan2(-a, n));
      TempDir d;
      UseGrid(d, {{"100", hfssfix::OnePointOneDirection(4.5e-3, 2.5e-5, P, T)},   // corner: fallback
                  {"200", OnePoint(0., 2.5e-5, AwayG(1e-3))},                      // normal inset, ray valid
                  {"400", OnePoint(0., 2.5e-5, TowardG(1e-3))}});                  // ray start
      World w(kExact); ScriptedEngine eng;
      BBSimOpBoundaryProcess::ResetRunTotals();
      CHECK(BBSimOpBoundaryProcess::RunTotals() == (Totals{0, 0, 0}));
      const double freqs[] = {100., 200., 400.};
      const Totals expect[] = {{1, 1, 0}, {2, 1, 1}, {3, 1, 1}};
      for (int i = 0; i < 3; ++i) {
        Script(eng, {0.0, 0.5, 0.5});
        CHECK(w.Shoot({-3 * mm, 0, 0}, X, Z, EnergyEV(freqs[i]), eng).status == S::kBBRDiffractionTransmit);
        CHECK(BBSimOpBoundaryProcess::RunTotals() == expect[i]);
      }
      BBSimOpBoundaryProcess::ResetRunTotals();
      CHECK(BBSimOpBoundaryProcess::RunTotals() == (Totals{0, 0, 0}));
    }},
    {"exit_start_summary_zero", [] {
      // The fallback count is printed when it is 0 too.
      TempDir d; UseData(d, hfssfix::WrapKeys());
      World w; ScriptedEngine eng;
      Capture cap;
      G4iosSetDestination(&cap);
      for (int i = 0; i < 100; ++i) {
        Script(eng, {0.0, 0.5, 0.5});
        w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      }
      G4iosSetDestination(nullptr);
      CHECK(cap.out.find("[BBR] diffraction events=100 T_obs=1 exit_start_fallback=0 exit_start_normal=0\n") !=
            std::string::npos);
    }},
    {"polarization_filter", [] {
      TempDir d; UseData(d, hfssfix::WrapKeys());
      World w; ScriptedEngine eng;
      const G4ThreeVector start(-3 * mm, 0, 0);
      // Normal incidence -> key (0,180), T0 = 1, T1 = 0; e_theta_in = -z, e_phi_in = -y.
      Script(eng, {0.999999, 0.5, 0.5});
      CHECK(w.Shoot(start, X, Z, E500, eng).status == S::kBBRDiffractionTransmit);   // E across the gap
      Script(eng, {1e-9});
      auto r = w.Shoot(start, X, Y, E500, eng);                                        // E along the long axis
      CHECK(r.status == S::kBBRDiffractionReflect);
      CHECK_VEC(r.dir, -X, 1e-12);
      // pol parallel to k: E_theta = E_phi = 1/sqrt2 -> T = (T0 + T1)/2 = 0.5
      Script(eng, {0.4999, 0.5, 0.5});
      CHECK(w.Shoot(start, X, X, E500, eng).status == S::kBBRDiffractionTransmit);
      Script(eng, {0.5001});
      CHECK(w.Shoot(start, X, X, E500, eng).status == S::kBBRDiffractionReflect);
      // 45 deg, k in the x-z plane -> key (0,135), T0 = 1, T1 = 0: e_theta_in lies in the
      // plane of incidence, e_phi_in = -y.
      const G4ThreeVector k(c45, 0, s45);
      Script(eng, {0.999999, 0.5, 0.5});
      CHECK(w.Shoot(kHitMinusX - 1 * mm * k, k, G4ThreeVector(c45, 0, -s45), E500, eng).status == S::kBBRDiffractionTransmit);
      Script(eng, {1e-9});
      CHECK(w.Shoot(kHitMinusX - 1 * mm * k, k, Y, E500, eng).status == S::kBBRDiffractionReflect);
    }},
    {"frequency_recorded", [] {
      TempDir d;
      hfssfix::WriteDataset(d.path(), "crack_100GHz", hfssfix::WrapKeys("100"));
      hfssfix::WriteDataset(d.path(), "crack_1000GHz", hfssfix::WrapKeys("1000"));
      BBRConfigManager::SetDataDir(d.path().string());
      World w; ScriptedEngine eng;
      Script(eng, {1e-9});
      CHECK_NEAR(w.Shoot({-3 * mm, 0, 0}, X, Y, EnergyEV(200.), eng).freqGHz, 100., 0.);   // log-nearest
      Script(eng, {1e-9});
      CHECK_NEAR(w.Shoot({-3 * mm, 0, 0}, X, Y, EnergyEV(400.), eng).freqGHz, 1000., 0.);  // not linear-nearest
    }},
    {"rng_budget", [] {
      TempDir d; UseData(d, hfssfix::WrapKeys());
      World w; ScriptedEngine eng;
      Script(eng, {0.999999, 0.5, 0.5});
      auto r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng);
      CHECK(r.status == S::kBBRDiffractionTransmit && r.uniformsUsed == 3 && eng.Remaining() == 0);
      Script(eng, {1e-9});
      r = w.Shoot({-3 * mm, 0, 0}, X, Y, E500, eng);
      CHECK(r.status == S::kBBRDiffractionReflect && r.uniformsUsed == 1 && eng.Remaining() == 0);
      // Not a geometry-boundary step: passed through to the stock process, no uniform, no lookup.
      Script(eng, {0.5});
      r = w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng, fAlongStepDoItProc);
      CHECK(r.status == S::kBBRNone && r.uniformsUsed == 0 && r.freqGHz == -1.);
      CHECK(w.Process().GetStatus() == NotAtBoundary);
    }},
    {"snap_and_rotation", [] {
      TempDir d; UseData(d, hfssfix::WrapKeys());
      ScriptedEngine eng;
      {
        World w;
        // k_y = +1e-13: kt = -1e-13 is snapped to +0 -> phi_raw = +180 -> sy = +1, sx = -1.
        const G4ThreeVector k1 = G4ThreeVector(c45, 1e-13, -s45).unit();
        Script(eng, {0.0, 0.5, 0.5});
        CHECK_VEC(w.Shoot(kHitMinusX - 1 * mm * k1, k1, InPlane(k1), E500, eng).pos,
                  ExpectedStart(G4ThreeVector(2 * mm, 1 * mm, -0.01 * mm), FFDir(30, 60, X, Y, -Z), X), 1e-9);
        // k_y = +1e-11 is not snapped -> phi_raw = -180 -> sy = -1, sx = -1.
        const G4ThreeVector k2 = G4ThreeVector(c45, 1e-11, -s45).unit();
        Script(eng, {0.0, 0.5, 0.5});
        CHECK_VEC(w.Shoot(kHitMinusX - 1 * mm * k2, k2, InPlane(k2), E500, eng).pos,
                  ExpectedStart(G4ThreeVector(2 * mm, -1 * mm, -0.01 * mm), FFDir(30, 60, X, -Y, -Z), X), 1e-9);
        // Entry through the +x face: the normal flips to -x, theta_hat/phi_hat do not.
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot({3 * mm, 0, 0}, -X, Z, E500, eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(-2 * mm, 1 * mm, 0.01 * mm), -X, -X), 1e-9);
        CHECK_VEC(r.dir, -X, 1e-12);
      }
      {
        // Active Rx(+90 deg): theta_hat (long) = world +z, phi_hat (gap) = world -y.
        G4RotationMatrix rx; rx.rotateX(90 * deg);
        World w({2 * mm, 5 * mm, 0.026 * mm}, rx);
        Script(eng, {0.999999, 0.5, 0.5});
        const auto r = w.Shoot({-3 * mm, 0, 0}, X, Y, E500, eng);   // world y is across the rotated gap
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(2 * mm, -0.01 * mm, 1 * mm), X, X), 1e-9);
        Script(eng, {1e-9});
        CHECK(w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng).status == S::kBBRDiffractionReflect);
      }
    }},
    {"fatal_bbr004_bbr005", [] {
      TempDir d; UseData(d, hfssfix::WrapKeys("500", "0.01"));   // exit point at Y = 10 mm
      ScriptedEngine eng;
      {
        World w;   // half-length 5 mm along y: the exit sample lies outside the crack
        Script(eng, {0.0, 0.5, 0.5});
        ExpectG4Exception("BBR005", [&] { w.Shoot({-3 * mm, 0, 0}, X, Z, E500, eng); }, "BBSimOpBoundaryProcess");
      }
      {
        World thin({1.5e-8 * mm, 5 * mm, 0.026 * mm});   // half-length 1.5e-8 mm <= inset 2e-8 mm (it passed at 10x)
        Script(eng, {0.0, 0.5, 0.5});
        ExpectG4Exception("BBR004", [&] { thin.Shoot({-1 * mm, 0, 0}, X, Z, E500, eng); }, "BBSimOpBoundaryProcess");
      }
    }},
    {"mirror_equivariance", [] {
      // The parallel plate is symmetric under y -> -y and z -> -z. With the same uniforms,
      // a mirrored input must give the same decision and the mirrored direction, exit
      // position (crack centred at the origin) and polarization (up to sign).
      TempDir d; UseData(d, hfssfix::NineKeyGrid());
      World w; ScriptedEngine eng;
      CLHEP::MixMaxRng gen(12345);
      auto My = [](G4ThreeVector v) { return G4ThreeVector(v.x(), -v.y(), v.z()); };
      auto Mz = [](G4ThreeVector v) { return G4ThreeVector(v.x(), v.y(), -v.z()); };
      const double Ud[] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};
      int nCase = 0, nBad = 0, nTransmit = 0;
      for (int i = 0; i < 400; ++i) {
        G4ThreeVector k;
        do { k = G4ThreeVector(0.2 + 0.8 * gen.flat(), 2 * gen.flat() - 1, 2 * gen.flat() - 1).unit(); }
        while (std::abs(k.y()) <= 1e-6 || std::abs(k.z()) <= 1e-6);
        const G4ThreeVector a = k.orthogonal().unit(), b = k.cross(a).unit();
        const double ph = CLHEP::twopi * gen.flat();
        const G4ThreeVector p = std::cos(ph) * a + std::sin(ph) * b;
        for (double ud : Ud) {
          const double u2 = gen.flat(), u3 = gen.flat();
          Script(eng, {ud, u2, u3});
          const auto o = w.Shoot(kHitMinusX - 1 * mm * k, k, p, E500, eng);
          nTransmit += (o.status == S::kBBRDiffractionTransmit);
          for (int m = 0; m < 3; ++m) {
            auto M = [&](G4ThreeVector v) { return m == 0 ? My(v) : (m == 1 ? Mz(v) : My(Mz(v))); };
            Script(eng, {ud, u2, u3});
            const auto r = w.Shoot(kHitMinusX - 1 * mm * M(k), M(k), M(p), E500, eng);
            ++nCase;
            bool ok = r.status == o.status && (r.dir - M(o.dir)).mag() < 1e-9 &&
                      std::abs(std::abs(r.pol.dot(M(o.pol))) - 1.) < 1e-9;
            if (o.status == S::kBBRDiffractionTransmit) ok = ok && (r.pos - M(o.pos)).mag() < 1e-9 * mm;
            if (!ok && nBad++ < 3)
              std::fprintf(stderr, "mismatch m=%d k=(%.4f,%.4f,%.4f) status %d/%d\n", m, k.x(), k.y(), k.z(),
                           int(o.status), int(r.status));
          }
        }
      }
      std::printf("cases=%d transmitted=%d of %d mismatches=%d\n", nCase, nTransmit, nCase / 3, nBad);
      CHECK(nCase == 8400);
      CHECK(nBad == 0);
      CHECK(nTransmit > 0 && nTransmit < nCase / 3);   // both branches exercised
    }},
    {"side_face_entry_documented", [] {
      // D10 characterization: entry through a side (+y) face is not detected; the photon
      // is treated as an axial entry and relocated to an x face. The geometry rule "only
      // the +-x faces of a vacuum_wg volume may be exposed" is in BBSimOpBoundaryProcess.hh.
      TempDir d; UseData(d, hfssfix::WrapKeys());
      World w; ScriptedEngine eng;
      for (double kx : {0.3, -0.3}) {
        const G4ThreeVector k = G4ThreeVector(kx, -0.95, 0).unit();
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot({0, 6 * mm, 0}, k, Z, E500, eng);
        CHECK_NEAR(r.hit.y(), 5 * mm, 1e-9);                      // entered through the +y face
        CHECK(r.status == S::kBBRDiffractionTransmit);            // key (90,135), T = 1
        const double sgn = kx > 0 ? 1. : -1.;                     // normal flipped along k
        CHECK_VEC(r.pos, ExpectedStart(G4ThreeVector(sgn * 2 * mm, 1 * mm, 0.01 * mm),
                                       FFDir(60, 90, sgn * X, Y, Z), sgn * X), 1e-9);
      }
    }},
    {"reflectance_decision", [] {
      // Flat R = 0.3: U <= R reflects, U > R absorbs (energy deposited, track killed); one uniform.
      World w({2 * mm, 5 * mm, 0.026 * mm}, {}, Flat30()); ScriptedEngine eng;
      const G4double E = 2.07e-3;
      Script(eng, {0.3});
      auto r = w.Shoot(kSlabStart, k30, Z, E, eng);
      CHECK(r.status == S::kBBRReflect && r.uniformsUsed == 1 && r.trackStatus == fAlive);
      Script(eng, {0.3000001});
      r = w.Shoot(kSlabStart, k30, Z, E, eng);
      CHECK(r.status == S::kBBRAbsorb && r.uniformsUsed == 1 && r.trackStatus == fStopAndKill);
      CHECK_REL(r.edep, E * eV, 1e-15);
      CHECK(Handler().CountWarnings("BBR006") == 0);   // the navigator gave a valid normal
    }},
    {"reflectance_specular", [] {
      World w({2 * mm, 5 * mm, 0.026 * mm}, {}, Flat30()); ScriptedEngine eng;
      const double c30 = std::cos(30 * deg), s30 = std::sin(30 * deg);
      Script(eng, {0.1});
      auto r = w.Shoot(kSlabStart, k30, Z, 2.07e-3, eng);            // s-polarization
      CHECK(r.status == S::kBBRReflect);
      CHECK_VEC(r.dir, G4ThreeVector(-c30, s30, 0), 1e-9);
      CHECK_VEC(r.pol, Z, 1e-12);
      Script(eng, {0.1});
      r = w.Shoot(kSlabStart, k30, G4ThreeVector(-s30, c30, 0), 2.07e-3, eng);   // p-polarization
      CHECK_AXIS(r.pol, G4ThreeVector(s30, c30, 0));
    }},
    {"reflectance_crack_wall", [] {
      // A photon inside the crack heading +z strikes the crack's +z wall and enters the slab
      // (its mother). It must reflect about the wall normal z, not about the slab face x
      // (the 2026-09-16 bug; caught only because Shoot drives a real navigated step).
      World w({2 * mm, 5 * mm, 0.026 * mm}, {}, Flat30()); ScriptedEngine eng;
      Script(eng, {0.1});
      const auto r = w.Shoot({0, 0, 0}, Z, X, 2.07e-3, eng);
      CHECK_NEAR(r.hit.z(), 0.026 * mm, 1e-9);
      CHECK(r.status == S::kBBRReflect);
      CHECK_VEC(r.dir, -Z, 1e-9);
      CHECK(Handler().CountWarnings("BBR006") == 0);
    }},
    {"reflectance_table_interp", [] {
      // Two-point table R(1e-3 eV) = 0, R(3e-3 eV) = 1: R(2e-3 eV) = 0.5 (linear in E).
      World w({2 * mm, 5 * mm, 0.026 * mm}, {}, Mirror("Rlin", {1e-3 * eV, 3e-3 * eV}, {0., 1.}));
      ScriptedEngine eng;
      Script(eng, {0.49});
      CHECK(w.Shoot(kSlabStart, k30, Z, 2e-3, eng).status == S::kBBRReflect);
      Script(eng, {0.51});
      CHECK(w.Shoot(kSlabStart, k30, Z, 2e-3, eng).status == S::kBBRAbsorb);
    }},
  });
}
