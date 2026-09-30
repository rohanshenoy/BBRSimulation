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
#include "CLHEP/Random/MixMaxRng.h"

#include <cmath>

using namespace bbrtest;
using wrapperworld::EnergyEV;
using wrapperworld::World;
using S = BBSimOpBoundaryProcess;

namespace {
const G4ThreeVector X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1);
const double c45 = std::cos(45 * deg), s45 = std::sin(45 * deg);
const G4ThreeVector kHitMinusX(-2 * mm, 0, 0);   // centre of the crack's -x face
const double E500 = EnergyEV(500.);

double Inset() { return 10. * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(); }

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

// Direction of far-field row (Phi, Theta) in the folded frame: sinT cosP n + sinT sinP tf + cosT pf.
G4ThreeVector FFDir(double P, double T, const G4ThreeVector& n, const G4ThreeVector& tf, const G4ThreeVector& pf) {
  P *= deg; T *= deg;
  return std::sin(T) * std::cos(P) * n + std::sin(T) * std::sin(P) * tf + std::cos(T) * pf;
}
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
        CHECK_VEC(r.pos, G4ThreeVector(2 * mm - Inset(), c.sy * 1 * mm, c.sx * 0.01 * mm), 1e-9);
        CHECK_VEC(r.dir, FFDir(c.P, c.T, X, c.sy * Y, c.sx * Z), 1e-9);
      }
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
                  G4ThreeVector(2 * mm - Inset(), 1 * mm, -0.01 * mm), 1e-9);
        // k_y = +1e-11 is not snapped -> phi_raw = -180 -> sy = -1, sx = -1.
        const G4ThreeVector k2 = G4ThreeVector(c45, 1e-11, -s45).unit();
        Script(eng, {0.0, 0.5, 0.5});
        CHECK_VEC(w.Shoot(kHitMinusX - 1 * mm * k2, k2, InPlane(k2), E500, eng).pos,
                  G4ThreeVector(2 * mm - Inset(), -1 * mm, -0.01 * mm), 1e-9);
        // Entry through the +x face: the normal flips to -x, theta_hat/phi_hat do not.
        Script(eng, {0.0, 0.5, 0.5});
        const auto r = w.Shoot({3 * mm, 0, 0}, -X, Z, E500, eng);
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.pos, G4ThreeVector(-2 * mm + Inset(), 1 * mm, 0.01 * mm), 1e-9);
        CHECK_VEC(r.dir, -X, 1e-12);
      }
      {
        // Active Rx(+90 deg): theta_hat (long) = world +z, phi_hat (gap) = world -y.
        G4RotationMatrix rx; rx.rotateX(90 * deg);
        World w({2 * mm, 5 * mm, 0.026 * mm}, rx);
        Script(eng, {0.999999, 0.5, 0.5});
        const auto r = w.Shoot({-3 * mm, 0, 0}, X, Y, E500, eng);   // world y is across the rotated gap
        CHECK(r.status == S::kBBRDiffractionTransmit);
        CHECK_VEC(r.pos, G4ThreeVector(2 * mm - Inset(), -0.01 * mm, 1 * mm), 1e-9);
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
        World thin({5e-9 * mm, 5 * mm, 0.026 * mm});   // half-length 5e-9 mm <= inset 1e-8 mm
        Script(eng, {0.0, 0.5, 0.5});
        ExpectG4Exception("BBR004", [&] { thin.Shoot({-1 * mm, 0, 0}, X, Z, E500, eng); }, "BBSimOpBoundaryProcess");
      }
    }},
    {"mirror_equivariance", [] {
      // The parallel plate is symmetric under y -> -y and z -> -z. With the same uniforms,
      // a mirrored input must give the same decision and the mirrored direction, exit
      // position (crack centred at the origin) and polarization (up to sign).
      TempDir d; UseData(d, hfssfix::NineKey());
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
        CHECK_VEC(r.pos, G4ThreeVector(sgn * (2 * mm - Inset()), 1 * mm, 0.01 * mm), 1e-9);
      }
    }},
  });
}
