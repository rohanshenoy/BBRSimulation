// testEmission.cc — GetBBSpecCDF (Planck photon-number CDF), ThermalSurface
// (exact events under CLHEP::NonRandomEngine, seeded statistics, guards
// BBR017-BBR019), GeometricSurface areas, and the band-loss fractions behind
// the examples' BBR021 warning.
#include "BBRTestSupport.hh"

#include "BBEvt.hh"
#include "GeometricSurface.hh"
#include "GetBBSpecCDF.hh"
#include "ThermalSurface.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/NonRandomEngine.h"

#include <cmath>
#include <vector>

namespace {

const double kEmin = 4.14e-5, kEmax = 8.27e-2;  // eV, the examples' band
const double kB    = 8.6173e-5;                 // eV/K, as GetBBSpecCDF.cc
const double kPi   = 3.14159265358979323846;
const double kUStar = 1.5936242600;             // root of 2(1 - e^-u) = u

void Seq(CLHEP::NonRandomEngine& eng, std::vector<double> v) {
  eng.setRandomSequence(v.data(), static_cast<int>(v.size()));
}

// The test-world emitter: 1 x 20 x 20 mm at (-50, 0, 0) mm, outward.
void TestWorldBox(ThermalSurface& ts, G4bool outward = true) {
  ts.BBSpecCDF.initialize(4., kEmin, kEmax);
  ts.AddBoxSurface(G4ThreeVector(-50. * mm, 0., 0.), 1. * mm, 20. * mm, 20. * mm, outward);
}

// Fraction of the photon-number spectrum u^2/(e^u - 1) outside
// [kEmin, kEmax] at temperature T: Simpson below u_min, the series
// sum_m e^{-m b}(b^2/m + 2b/m^2 + 2/m^3) above u_max, over 2 zeta(3).
double BandLoss(double T) {
  const double a = kEmin / (kB * T), b = kEmax / (kB * T);
  auto f = [](double u) { return u > 0. ? u * u / std::expm1(u) : 0.; };
  const int n = 20000;
  const double h = a / n;
  double s = f(0.) + f(a);
  for (int i = 1; i < n; ++i) s += f(i * h) * (i % 2 ? 4. : 2.);
  const double lower = s * h / 3.;
  double upper = 0.;
  for (int m = 1; m <= 60; ++m)
    upper += std::exp(-m * b) * (b * b / m + 2. * b / (m * m) + 2. / (double(m) * m * m));
  return (lower + upper) / (2. * 1.2020569031595942);
}

}  // namespace

int main(int argc, char** argv) {
  using namespace bbrtest;
  return RunCase(argc, argv, {
    {"planck_cdf_shape", [] {
      const double Ts[] = {4., 10., 20.};
      const double meanE[] = {9.3368612e-4, 2.3287732e-3, 4.6559233e-3};
      for (int t = 0; t < 3; ++t) {
        GetBBSpecCDF c;
        c.initialize(Ts[t], kEmin, kEmax);
        const std::size_t N = c.x.size();
        CHECK(N == 100001 && c.pdf.size() == N && c.cdf.size() == N);
        CHECK(c.cdf[0] == 0.);
        CHECK(c.cdf[N - 1] == 1.);
        bool mono = true;
        double norm = 0., mean = 0.;
        std::size_t imax = 0;
        for (std::size_t i = 1; i < N; ++i) {
          mono = mono && c.cdf[i] >= c.cdf[i - 1];
          const double dx = c.x[i] - c.x[i - 1];
          norm += 0.5 * (c.pdf[i] + c.pdf[i - 1]) * dx;
          mean += 0.5 * (c.x[i] * c.pdf[i] + c.x[i - 1] * c.pdf[i - 1]) * dx;
          if (c.pdf[i] > c.pdf[imax]) imax = i;
        }
        CHECK(mono);
        CHECK_NEAR(norm, 1., 1e-12);
        CHECK_NEAR(c.x[imax], kUStar * kB * Ts[t], c.x[1] - c.x[0]);
        CHECK_REL(mean, meanE[t], 1e-5);
      }
    }},
    {"planck_cdf_invalid", [] {  // D8
      ExpectG4Exception("BBR017", [] { GetBBSpecCDF c; c.initialize(0., kEmin, kEmax); }, "GetBBSpecCDF");
      ExpectG4Exception("BBR017", [] { GetBBSpecCDF c; c.initialize(4., 0., kEmax); }, "GetBBSpecCDF");
      ExpectG4Exception("BBR017", [] { GetBBSpecCDF c; c.initialize(4., kEmax, kEmin); }, "GetBBSpecCDF");
      ExpectG4Exception("BBR017", [] { GetBBSpecCDF c; c.initialize(std::nan(""), kEmin, kEmax); }, "GetBBSpecCDF");
    }},
    {"planck_band_loss", [] {
      // P3: the band is fixed at 10 GHz-20 THz. The examples warn (BBR021)
      // outside 2.2-117 K, where more than 1 % of the photon-number spectrum
      // falls outside it; these are the numbers behind those two literals.
      CHECK_REL(BandLoss(0.5), 0.1378267, 1e-4);
      CHECK_REL(BandLoss(4.), 2.881908e-3, 1e-4);
      CHECK_REL(BandLoss(150.), 3.862435e-2, 1e-4);
      CHECK(BandLoss(2.2) <= 0.01 && BandLoss(2.1) > 0.01);
      CHECK(BandLoss(117.) <= 0.01 && BandLoss(118.) > 0.01);
    }},
    {"genevt_nonrandom_exact", [] {
      // Draw order: uE, uSurface, uFaceSign, uTheta, uPhi, uFaceSelect, coordA, coordB.
      auto* eng = new CLHEP::NonRandomEngine;
      G4Random::setTheEngine(eng);
      ThermalSurface ts;
      TestWorldBox(ts);
      Seq(*eng, {0.5, 0.3, 0.7, 0.5, 0.25, 0.1, 0.75, 0.25});
      BBEvt e = ts.GenEvt();
      CHECK_REL(e.energy, 8.144136472552e-4, 1e-9);
      CHECK_VEC(e.position, G4ThreeVector(-49.5, 5., -5.) * mm, 1e-12);
      CHECK_VEC(e.direction, G4ThreeVector(1., 0., 1.) / std::sqrt(2.), 1e-12);
      eng->setNextRandom(0.123);  // unread if GenEvt took exactly 8 uniforms
      CHECK(eng->flat() == 0.123);
      // uFaceSign = 0.2 and theta = 0: the -x face, emitting along -x.
      Seq(*eng, {0.5, 0.3, 0.2, 0.0, 0.0, 0.1, 0.5, 0.5});
      e = ts.GenEvt();
      CHECK_VEC(e.position, G4ThreeVector(-50.5, 0., 0.) * mm, 1e-12);
      CHECK_VEC(e.direction, G4ThreeVector(-1., 0., 0.), 1e-12);
      // Face thresholds at 800/880 and 840/880: 0.93 -> a y face, 0.99 -> a z face.
      Seq(*eng, {0.5, 0.3, 0.7, 0.5, 0.0, 0.93, 0.75, 0.25});
      e = ts.GenEvt();
      CHECK_NEAR(std::abs(e.position.y()), 10. * mm, 1e-12);
      Seq(*eng, {0.5, 0.3, 0.2, 0.5, 0.0, 0.99, 0.75, 0.25});
      e = ts.GenEvt();
      CHECK_NEAR(std::abs(e.position.z()), 10. * mm, 1e-12);
      // Energy quantiles against the truncated Planck photon-number CDF (scipy).
      const double uE[] = {0.1, 0.5, 0.9};
      const double qE[] = {2.7806311e-4, 8.1441348e-4, 1.7447505e-3};
      for (int i = 0; i < 3; ++i) {
        Seq(*eng, {uE[i], 0.3, 0.7, 0.5, 0.0, 0.1, 0.5, 0.5});
        CHECK_REL(ts.GenEvt().energy, qE[i], 1e-5);
      }
      // in_out = false: same point, the normal component of the direction flips.
      ThermalSurface in;
      TestWorldBox(in, false);
      Seq(*eng, {0.5, 0.3, 0.7, 0.5, 0.25, 0.1, 0.75, 0.25});
      e = in.GenEvt();
      CHECK_VEC(e.position, G4ThreeVector(-49.5, 5., -5.) * mm, 1e-12);
      CHECK_VEC(e.direction, G4ThreeVector(-1., 0., 1.) / std::sqrt(2.), 1e-12);
    }},
    {"genevt_rotation", [] {  // D16: fixed-axis Z, then Y, then X
      auto* eng = new CLHEP::NonRandomEngine;
      G4Random::setTheEngine(eng);
      // +x face, theta = 0, local point (0.5, 5, -5) mm.
      const std::vector<double> s = {0.5, 0.3, 0.7, 0.0, 0.0, 0.1, 0.75, 0.25};
      ThermalSurface r1;
      r1.BBSpecCDF.initialize(4., kEmin, kEmax);
      r1.AddBoxSurface(G4ThreeVector(), 1. * mm, 20. * mm, 20. * mm, true, 90. * deg, 0., 0.);
      Seq(*eng, s);
      BBEvt e = r1.GenEvt();
      CHECK_VEC(e.position, G4ThreeVector(-5., 0.5, -5.) * mm, 1e-12);
      CHECK_VEC(e.direction, G4ThreeVector(0., 1., 0.), 1e-12);
      ThermalSurface r2;
      r2.BBSpecCDF.initialize(4., kEmin, kEmax);
      r2.AddBoxSurface(G4ThreeVector(), 1. * mm, 20. * mm, 20. * mm, true, 90. * deg, 90. * deg, 0.);
      Seq(*eng, s);
      e = r2.GenEvt();
      CHECK_VEC(e.position, G4ThreeVector(-5., 0.5, 5.) * mm, 1e-12);
      CHECK_VEC(e.direction, G4ThreeVector(0., 1., 0.), 1e-12);
    }},
    {"genevt_statistics", [] {
      // Uniform-in-theta emission (Chang's convention): <theta> = pi/4 and
      // <cos theta> = 2/pi, far from Lambert's 2/3. Tolerances, not pinned
      // MixMax numbers.
      G4Random::setTheEngine(new CLHEP::MixMaxRng);
      G4Random::setTheSeed(12345);
      ThermalSurface ts;
      TestWorldBox(ts);
      const G4ThreeVector C(-50. * mm, 0., 0.);
      const int N = 200000;
      long nface[3] = {0, 0, 0}, off = 0;
      double sumTh = 0., sumCos = 0.;
      for (int i = 0; i < N; ++i) {
        const BBEvt e = ts.GenEvt();
        const G4ThreeVector p = e.position - C;
        G4ThreeVector nrm;
        if (std::abs(std::abs(p.x()) - 0.5 * mm) < 1e-9) { ++nface[0]; nrm = G4ThreeVector(p.x() > 0 ? 1 : -1, 0, 0); }
        else if (std::abs(std::abs(p.y()) - 10. * mm) < 1e-9) { ++nface[1]; nrm = G4ThreeVector(0, p.y() > 0 ? 1 : -1, 0); }
        else if (std::abs(std::abs(p.z()) - 10. * mm) < 1e-9) { ++nface[2]; nrm = G4ThreeVector(0, 0, p.z() > 0 ? 1 : -1); }
        else { ++off; continue; }
        const double c = std::min(1., e.direction.unit().dot(nrm));
        sumCos += c;
        sumTh += std::acos(c);
      }
      CHECK(off == 0);
      const double sTh = std::sqrt(kPi * kPi / 48. / N);                 // ~1.0e-3
      const double sCos = std::sqrt((0.5 - 4. / (kPi * kPi)) / N);       // ~6.9e-4
      CHECK_NEAR(sumTh / N, kPi / 4., 4. * sTh);
      CHECK_NEAR(sumCos / N, 2. / kPi, 4. * sCos);
      CHECK(std::abs(sumCos / N - 2. / 3.) > 20. * sCos);
      const double pf[3] = {10. / 11., 1. / 22., 1. / 22.};
      for (int f = 0; f < 3; ++f)
        CHECK_NEAR(double(nface[f]) / N, pf[f], 4. * std::sqrt(pf[f] * (1. - pf[f]) / N));
      // Emissivity weighting: two equal boxes, emissivity 1 and 0.25.
      ThermalSurface two;
      two.BBSpecCDF.initialize(4., kEmin, kEmax);
      two.AddBoxSurface(G4ThreeVector(-100. * mm, 0., 0.), 10. * mm, 10. * mm, 10. * mm, true, 0., 0., 0., 1.0);
      two.AddBoxSurface(G4ThreeVector(+100. * mm, 0., 0.), 10. * mm, 10. * mm, 10. * mm, true, 0., 0., 0., 0.25);
      long nA = 0;
      for (int i = 0; i < N; ++i) nA += (two.GenEvt().position.x() < 0.);
      CHECK_NEAR(double(nA) / N, 0.80, 0.004);
      // Same seed, bit-identical first event.
      G4Random::setTheSeed(12345);
      const BBEvt e1 = ts.GenEvt();
      G4Random::setTheSeed(12345);
      const BBEvt e2 = ts.GenEvt();
      CHECK(e1.energy == e2.energy && e1.position == e2.position && e1.direction == e2.direction);
    }},
    {"genevt_guards", [] {  // D5, D12
      ExpectG4Exception("BBR018", [] {
        ThermalSurface ts;  // CDF never initialized (used to read pdf[-1])
        ts.AddBoxSurface(G4ThreeVector(), 1. * mm, 1. * mm, 1. * mm, true);
        ts.GenEvt();
      }, "ThermalSurface::GenEvt");
      ExpectG4Exception("BBR019", [] {
        ThermalSurface ts;
        ts.BBSpecCDF.initialize(4., kEmin, kEmax);
        ts.GenEvt();
      }, "ThermalSurface::GenEvt");
    }},
    {"geometric_surface_area", [] {
      GeometricSurface t;
      t.type = 1; t.radius = 1.; t.height = 1.;
      CHECK_REL(t.CalculateArea(), 4. * kPi, 1e-14);
      t.lid = true;
      CHECK_REL(t.CalculateArea(), 6. * kPi, 1e-14);
      GeometricSurface d;
      d.type = 2; d.radius = 2.; d.r_in = 1.;
      CHECK_REL(d.CalculateArea(), 3. * kPi, 1e-14);
      d.both_side = true;
      CHECK_REL(d.CalculateArea(), 6. * kPi, 1e-14);
      GeometricSurface s;
      s.type = 4; s.radius = 1.; s.theta1 = 0.; s.theta2 = kPi; s.phi1 = 0.; s.phi2 = 2. * kPi;
      CHECK_REL(s.CalculateArea(), 4. * kPi, 1e-14);
      GeometricSurface b;
      b.type = 3; b.Wx = 1.; b.Wy = 20.; b.Wz = 20.;
      CHECK(b.CalculateArea() == 880.);
      GeometricSurface u;
      u.type = 7; u.area = 5.;
      CHECK(u.CalculateArea() == 0. && u.area == 0.);
    }},
  });
}
