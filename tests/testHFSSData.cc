#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "BBRHFSSData.hh"
#include "BBRDatasetSidecar.hh"
#include "G4SystemOfUnits.hh"
#include "G4coutDestination.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "CLHEP/Random/MixMaxRng.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

using namespace bbrtest;
using hfssfix::Dataset;

namespace {
const G4ThreeVector N(1, 0, 0), TH(0, 1, 0), PH(0, 0, 1);   // normal, theta_hat, phi_hat
const double R2 = 1. / std::sqrt(2.);

// Collects what G4cout receives while installed (the loader logs the normalization there).
struct Capture : public G4coutDestination {
  std::string out;
  G4int ReceiveG4cout(const G4String& s) override { out += s; return 0; }
};

// Loads a dataset written into a fresh TempDir; the TempDir must outlive the object.
BBRHFSSData Load(TempDir& d, const Dataset& ds, const std::string& stem = "c_500GHz", double f = 500.) {
  const auto wg = hfssfix::WriteDataset(d.path(), stem, ds);
  return BBRHFSSData(wg.string(), stem, f);
}

// The dataset under <TempDir>/waveguides/c_500GHz with the sidecar built from p
// (default: SidecarFrom), loaded through BBRHFSSData with that sidecar (C1-C5).
BBRHFSSData LoadChecked(TempDir& d, const Dataset& ds, const hfssfix::SidecarParams* p = nullptr) {
  const std::string stem = "c_500GHz";
  const auto wg = hfssfix::WriteDataset(d.path(), stem, ds, false);
  const auto params = p ? *p : hfssfix::SidecarFrom(stem, ds);
  const auto sc = BBRDatasetSidecar::Parse(hfssfix::SidecarJson(params), "test", "c", stem, 500.);
  return BBRHFSSData(wg.string(), stem, 500., &sc);
}

// Keys with their own far-field rows ("Phi,Theta") and the same exit rows ("Y,Z") in both
// polarizations: T0 = 3/4, T1 = 1/4, rEtheta = 1 (Ephi=0) and rEphi = 1 (Ephi=1) at every far-field
// row, Ey (Ephi=0) and Ez (Ephi=1) at every exit point (the CDX-022 cases).
Dataset Keys(const std::vector<std::pair<std::string, std::vector<std::string>>>& ffByKey,
             const std::vector<std::string>& ex) {
  Dataset d{hfssfix::FF_HDR, hfssfix::WG_HDR, hfssfix::FF_HDR, hfssfix::WG_HDR};
  for (const auto& [k, ff] : ffByKey) {
    for (const auto& r : ff) {
      d.ff0 += "500GHz,0," + k + "," + r + ",0,0,1,0\n";
      d.ff1 += "500GHz,1," + k + "," + r + ",1,0,0,0\n";
    }
    for (const auto& r : ex) {
      d.wg0 += "500GHz,0," + k + ",3,4,0," + r + ",0,1,0,0,0,0\n";
      d.wg1 += "500GHz,1," + k + ",1,4,0," + r + ",0,0,1,0,0,0\n";
    }
  }
  return d;
}
Dataset Rows(const std::vector<std::string>& ff, const std::vector<std::string>& ex) { return Keys({{"0,180", ff}}, ex); }

// text with every `from` replaced by `to`; throws if `from` is absent, so a stale edit cannot pass.
std::string ReplaceAll(std::string text, const std::string& from, const std::string& to) {
  if (text.find(from) == std::string::npos) throw std::runtime_error("ReplaceAll: '" + from + "' not in the text");
  for (auto p = text.find(from); p != std::string::npos; p = text.find(from, p + to.size())) text.replace(p, from.size(), to);
  return text;
}

// A complete exit lattice at one key (0,180) with one far-field row and the same field at every exit
// point (Ey for Ephi=0, Ez for Ephi=1): the exit CDF is uniform, so row j of M is drawn with
// U = (j + 0.5) / M (the OR-020 clamp tests).
using YZ = std::pair<double, double>;
Dataset UniformLattice(const std::vector<YZ>& pts) {
  const auto num = [](double v) { std::ostringstream o; o << std::setprecision(17) << v; return o.str(); };
  std::vector<std::string> ex;
  for (const auto& [y, z] : pts) ex.push_back(num(y) + "," + num(z));
  return Rows({"0,90"}, ex);
}
// Where SampleExitPosition puts row j of M, with the exit face's centre at 0 and the axes (Y, Z).
G4ThreeVector DrawRow(const BBRHFSSData& h, ScriptedEngine& eng, std::size_t j, std::size_t M) {
  eng.Push((double(j) + 0.5) / double(M));
  return h.SampleExitPosition(1, 0, 0, 180, G4ThreeVector(), TH, PH);
}
std::size_t Count(const std::string& text, const std::string& what) {
  std::size_t n = 0;
  for (auto q = text.find(what); q != std::string::npos; q = text.find(what, q + what.size())) ++n;
  return n;
}
std::string ClampLine(const std::string& stem, int n) {
  return "[BBR] HFSS " + stem + ": " + std::to_string(n) + " exit point(s) clamped 1 nm inside the declared cross-section\n";
}

// --- schema 2 (hfssfix::NineKeyGrid2) ----------------------------------------------
const std::string kS2 = "c_500GHz";
const double kPi = 3.14159265358979323846;

// The dataset written under <TempDir>/waveguides/c_500GHz with the 2.0 sidecar of p, loaded with it.
BBRHFSSData LoadV2(TempDir& d, const Dataset& ds, const hfssfix::Schema2Params& p) {
  const auto wg = hfssfix::WriteDataset2(d.path(), kS2, ds, p);
  const auto sc = BBRDatasetSidecar::Parse(hfssfix::SidecarJson2(p), "test", "c", kS2, 500.);
  return BBRHFSSData(wg.string(), kS2, 500., &sc);
}

// csv with f(fields, row index within the key) applied to every data row of key (phi, theta).
std::string EditKey(const std::string& csv, double phi, double theta,
                    const std::function<void(std::vector<std::string>&, int)>& f) {
  std::istringstream in(csv);
  std::string line, out;
  bool header = true;
  int n = 0;
  while (std::getline(in, line)) {
    if (!header && !line.empty()) {
      std::vector<std::string> v;
      std::stringstream ss(line);
      std::string t;
      while (std::getline(ss, t, ',')) v.push_back(t);
      if (std::stod(v[2]) == phi && std::stod(v[3]) == theta) {
        f(v, n++);
        line.clear();
        for (std::size_t i = 0; i < v.size(); ++i) line += (i ? "," : "") + v[i];
      }
    }
    header = false;
    out += line + "\n";
  }
  if (n == 0) throw std::runtime_error("EditKey: no row of the key");
  return out;
}
std::string R17(double v) { return hfssfix::detail::Fmt(v); }

// NineKeyGrid2's stored field value of lattice point i (i = 3 iy + iz, y {-1, 0, 1} mm outer,
// z {-25, 0, 25} um inner) in polarization e, slot f = 0..5 (CSV positions 9..14).
double F2(int e, int i, int f) { return (e == 0 ? i + 1 + f : i + 6 - f) * 1e-3; }
// The 2.0 positions: Ex = (9, 10), Ey = (11, 12) as (re, im).
std::complex<double> Ex2(int e, int i) { return {F2(e, i, 0), F2(e, i, 1)}; }
std::complex<double> Ey2(int e, int i) { return {F2(e, i, 2), F2(e, i, 3)}; }
// The trapezoid weight of point i on the 3 x 3 lattice: 1/2 per wall row by index, a corner 1/4.
double W2(int i) { return (i / 3 == 1 ? 1. : 0.5) * (i % 3 == 1 ? 1. : 0.5); }

using P9 = std::array<double, 9>;
P9 Normalized(P9 p) {
  double s = 0.;
  for (double v : p) s += v;
  for (double& v : p) v /= s;
  return p;
}
// Predicted exit-point fractions for the mixture (a, b) of the stored fields: w_i |a E_t0 + b E_t1|^2.
P9 ExitPrediction(double a, double b, bool weighted) {
  P9 p{};
  for (int i = 0; i < 9; ++i)
    p[i] = (weighted ? W2(i) : 1.) * (std::norm(a * Ex2(0, i) + b * Ex2(1, i)) + std::norm(a * Ey2(0, i) + b * Ey2(1, i)));
  return Normalized(p);
}

// n exit positions drawn at (phi, theta) for (Et, Ep) with a seeded MixMax engine, counted per
// lattice point (the drawn position rounded onto the 1 mm x 25 um lattice, so a point clamped
// 1 nm inside the wall counts as its node). maxAbsZ_mm gets the largest |z| drawn.
std::array<long, 9> DrawExits(const BBRHFSSData& h, double Et, double Ep, double phi, double theta, long n,
                              long seed, double* maxAbsZ_mm = nullptr) {
  G4Random::setTheEngine(new CLHEP::MixMaxRng(seed));
  std::array<long, 9> c{};
  double zmax = 0.;
  for (long k = 0; k < n; ++k) {
    const auto r = h.SampleExitPosition(Et, Ep, phi, theta, G4ThreeVector(), TH, PH);
    const long iy = std::lround(r.y() / CLHEP::mm) + 1, iz = std::lround(r.z() / (0.025 * CLHEP::mm)) + 1;
    if (iy >= 0 && iy < 3 && iz >= 0 && iz < 3) ++c[3 * iy + iz];
    zmax = std::max(zmax, std::abs(r.z() / CLHEP::mm));
  }
  if (maxAbsZ_mm) *maxAbsZ_mm = zmax;
  return c;
}
// (observed fraction - p) / sigma of a binomial count.
double Zscore(long count, long n, double p) { return (double(count) / n - p) / std::sqrt(p * (1. - p) / n); }
// Checks every point's count within 4 sigma of p; prints the counts; returns the largest |z|.
double CheckCounts(const char* what, const std::array<long, 9>& c, long n, const P9& p) {
  long sum = 0;
  double zmax = 0.;
  for (int i = 0; i < 9; ++i) {
    const double z = Zscore(c[i], n, p[i]);
    std::printf("%s point %d: %ld of %ld, predicted %.6f, z = %+.3f\n", what, i, c[i], n, p[i], z);
    CHECK(std::abs(z) <= 4.);
    zmax = std::max(zmax, std::abs(z));
    sum += c[i];
  }
  CHECK(sum == n);   // every draw landed on a lattice point
  return zmax;
}

// NineKeyGrid's far-field row amplitudes (rEphi_re, rEphi_im, rEtheta_re, rEtheta_im) per polarization,
// on Phi {-40, 10, 60} x Theta {40, 75, 110}; row r = 3 itheta + iphi.
const double kFFPhi[3] = {-40., 10., 60.}, kFFTheta[3] = {40., 75., 110.};
std::array<double, 4> FF(int e, int r) {
  const double phi = kFFPhi[r % 3], theta = kFFTheta[r / 3];
  if (phi == 10. && theta == 40.) return e ? std::array<double, 4>{0.3, 0.1, 0.8, -0.2} : std::array<double, 4>{0.9, 0.0, 0.1, 0.4};
  if (phi == -40. && theta == 75.) return e ? std::array<double, 4>{-0.5, 0.2, 0.1, 0.6} : std::array<double, 4>{0.2, -0.3, 0.7, 0.1};
  if (phi == 60. && theta == 110.) return e ? std::array<double, 4>{0.4, 0.4, -0.3, 0.2} : std::array<double, 4>{0.5, 0.1, -0.6, 0.3};
  return e ? std::array<double, 4>{0.02, -0.03, 0.04, 0.01} : std::array<double, 4>{0.05, 0.02, 0.03, -0.01};
}
P9 DirectionPrediction(double a, double b) {
  P9 p{};
  for (int r = 0; r < 9; ++r) {
    const auto f0 = FF(0, r), f1 = FF(1, r);
    double s = 0.;
    for (int k = 0; k < 4; ++k) s += (a * f0[k] + b * f1[k]) * (a * f0[k] + b * f1[k]);
    p[r] = std::sin(kFFTheta[r / 3] * kPi / 180.) * s;
  }
  return Normalized(p);
}
std::array<long, 9> DrawDirections(const BBRHFSSData& h, double Et, double Ep, double phi, double theta, long n, long seed) {
  G4Random::setTheEngine(new CLHEP::MixMaxRng(seed));
  std::array<long, 9> c{};
  G4ThreeVector pol;
  for (long k = 0; k < n; ++k) {
    const auto dir = h.SampleOutgoingDirection(Et, Ep, phi, theta, PH, TH, N, pol);
    const double T = std::acos(std::max(-1., std::min(1., dir.z()))) * 180. / kPi;
    const double P = std::atan2(dir.y(), dir.x()) * 180. / kPi;
    for (int r = 0; r < 9; ++r)
      if (std::abs(T - kFFTheta[r / 3]) < 1e-6 && std::abs(P - kFFPhi[r % 3]) < 1e-6) ++c[r];
  }
  return c;
}

}  // namespace

int main(int argc, char** argv) {
  return RunCase(argc, argv, {
    {"parse_frequency", [] {
      using H = BBRHFSSData;
      CHECK_NEAR(H::ParseFrequencyGHz("500GHz"), 500., 0.);
      CHECK_NEAR(H::ParseFrequencyGHz("0.5THz"), 500., 1e-12);
      CHECK_NEAR(H::ParseFrequencyGHz("500000MHz"), 500., 1e-12);
      CHECK_NEAR(H::ParseFrequencyGHz(" 500 ghz "), 500., 0.);
      CHECK_NEAR(H::ParseFrequencyGHz("1.5e3GHz"), 1500., 0.);
      CHECK(H::ParseFrequencyGHz("500") < 0);
      CHECK(H::ParseFrequencyGHz("500Hz") < 0);
      CHECK(H::ParseFrequencyGHz("500 kHz") < 0);
      CHECK(H::ParseFrequencyGHz("-5GHz") < 0);
      CHECK(H::ParseFrequencyGHz("0GHz") < 0);
      CHECK(H::ParseFrequencyGHz("nanGHz") < 0);
      CHECK(H::ParseFrequencyGHz("") < 0);
      CHECK(H::ParseFrequencyGHz("GHz") < 0);
    }},
    {"mini_loads", [] {
      TempDir d; auto h = Load(d, hfssfix::Mini500());
      CHECK_NEAR(h.GetFrequencyGHz(), 500., 0.);
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 0.75, 1e-15);
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 180), 0.25, 1e-15);
    }},
    {"transmittance_mixing", [] {
      // T0 = 0.8, T1 = 0.2 at (0,180); orthogonal exit fields, so no cross term.
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,0.8,1,0,0,0,0,1,0,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.2,1,0,0,0,0,0,1,0,0,0\n";
      TempDir d; auto h = Load(d, ds);
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 0.8, 1e-15);
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 180), 0.2, 1e-15);
      CHECK_NEAR(h.GetTransmittance(R2, R2, 0, 180), 0.5, 1e-15);
      CHECK_NEAR(h.GetTransmittance(0.6, 0.8, 0, 180), 0.36 * 0.8 + 0.64 * 0.2, 1e-15);
      CHECK_NEAR(h.GetTransmittance(2, 0, 0, 180), 1.0, 0.);   // clamp
    }},
    {"load_cap_first_row", [] {
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,1.2,1,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,0.123,1,0,0.001,0,0,1,0,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.5,0,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,0.5,0,0,0.001,0,0,0,1,0,0,0\n";
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir d; auto h = Load(d, ds);
      G4iosSetDestination(nullptr);
      CHECK(cap.out.find("normalized to 1") != std::string::npos);   // the normalization is logged
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 1.0, 0.);    // raw 1.2 divided by lambda_max = 1.2
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 180), 0.0, 0.);    // IngoingPower 0 -> 0 (0 / 1.2)
    }},
    {"nearest_key", [] {
      TempDir d; auto h = Load(d, hfssfix::ThreeKey());
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 170), 0.8, 1e-15);   // nearest (0,180)
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 150), 1.0, 1e-15);   // nearest (0,135)
      CHECK_NEAR(h.GetTransmittance(1, 0, 46, 180), 0.3, 1e-15);  // L2 distance to (0,180) is 46, to (90,180) is 44, so (90,180) wins
      CHECK_NEAR(h.GetTransmittance(1, 0, 45, 180), 0.8, 1e-15);  // exact tie (0,180)/(90,180): std::map order keeps (0,180), as hfss.nearest_key
    }},
    {"direction_cdf", [] {
      // Three rows at Theta = 90 (sinT = 1), Phi = 0 / 30 / 60; Ephi=0 amplitudes
      // rEtheta = 1, 1 and (rEphi, rEtheta) = (1, 1) -> weights exactly 1, 1, 2 ->
      // CDF 0.25, 0.5, 1.0 (a rounded sqrt(2) literal would move the 0.25 edge
      // by ~5e-13, past the "just above" draw).
      Dataset ds = hfssfix::Mini500();
      ds.ff0 = std::string(hfssfix::FF_HDR) + "500GHz,0,0,180,0,90,0,0,1,0\n500GHz,0,0,180,30,90,0,0,1,0\n500GHz,0,0,180,60,90,1,0,1,0\n";
      // Ephi=1 (T1 = 1/4) radiates too, as every transmitting polarization must (BBR027); E_phi = 0 never sees it.
      ds.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90,1,0,0,0\n500GHz,1,0,180,30,90,0,0,0,0\n500GHz,1,0,180,60,90,0,0,0,0\n";
      TempDir d; auto h = Load(d, ds);
      ScriptedEngine eng{0.25, 0.2500000000001, 0.5, 0.75, 1.0};
      G4Random::setTheEngine(&eng);
      G4ThreeVector pol;
      const double c30 = std::cos(M_PI / 6), s30 = 0.5, c60 = 0.5, s60 = std::sin(M_PI / 3);
      CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), N, 1e-12);                       // U = 0.25 -> row 0
      CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), c30 * N + s30 * TH, 1e-9);      // just above -> row 1
      CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), c30 * N + s30 * TH, 1e-9);      // U = 0.5 -> row 1
      CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), c60 * N + s60 * TH, 1e-9);      // U = 0.75 -> row 2
      CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), c60 * N + s60 * TH, 1e-9);      // U = 1 -> row 2
      CHECK(eng.Draws() == 5);   // exactly one uniform per call
    }},
    {"zero_weight_rows", [] {
      // Row A in the exit-face plane (Phi = -90), row B at Theta = 0, row C straight out.
      Dataset ds = hfssfix::Mini500();
      ds.ff0 = std::string(hfssfix::FF_HDR) + "500GHz,0,0,180,-90,90,0,0,1,0\n500GHz,0,0,180,0,0,0,0,1,0\n500GHz,0,0,180,0,90,0,0,1,0\n";
      ds.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,-90,90,0,0,0,0\n500GHz,1,0,180,0,0,0,0,0,0\n500GHz,1,0,180,0,90,1,0,0,0\n";
      TempDir d; auto h = Load(d, ds);
      ScriptedEngine eng{1e-12, 0.5, 1.0};
      G4Random::setTheEngine(&eng);
      G4ThreeVector pol;
      for (int i = 0; i < 3; ++i)
        CHECK_VEC(h.SampleOutgoingDirection(1, 0, 0, 180, PH, TH, N, pol), N, 1e-12);   // always row C
    }},
    {"coherent_cross_term_in_direction", [] {
      // Row B: F0 = +1, F1 = -1 (theta components). With Et = Ep = 1/sqrt2 its weight is 0.
      // Row C: F0 = 2i along phi, F1 = 1 along theta, so its amplitude is (a, b) =
      // (1, 2i)/sqrt2 in (theta_hat, phi_hat), an ellipse with its major axis along
      // phi_hat; the real parts alone would point along theta_hat.
      Dataset ds = hfssfix::Mini500();
      ds.ff0 = std::string(hfssfix::FF_HDR) + "500GHz,0,0,180,0,90,0,0,1,0\n500GHz,0,0,180,60,90,0,2,0,0\n";
      ds.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90,0,0,-1,0\n500GHz,1,0,180,60,90,0,0,1,0\n";
      TempDir d; auto h = Load(d, ds);
      ScriptedEngine eng{0.01, 0.99};
      G4Random::setTheEngine(&eng);
      G4ThreeVector pol;
      const G4ThreeVector rowC = 0.5 * N + std::sin(M_PI / 3) * TH;
      const G4ThreeVector phiC = -std::sin(M_PI / 3) * N + 0.5 * TH;   // phi_hat at row C
      CHECK_VEC(h.SampleOutgoingDirection(R2, R2, 0, 180, PH, TH, N, pol), rowC, 1e-9);
      CHECK_AXIS(pol, phiC);                                            // the major axis
      CHECK_VEC(h.SampleOutgoingDirection(R2, R2, 0, 180, PH, TH, N, pol), rowC, 1e-9);
    }},
    {"exit_position", [] {
      // Exit points at Y = 1 mm and Y = -2 mm, Z = 10 um (X = 5 mm, which the sampler
      // ignores); fields (0,1,0) and (1,1,1) -> |E|^2 weights exactly 1 and 3 -> CDF 0.25, 1.
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,1,1,0.005,0.001,0,0,1,0,0,0,0\n500GHz,0,0,180,1,1,0.005,-0.002,0.00001,1,1,1,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,1,0.005,0.001,0,0,0,1,0,0,0\n500GHz,1,0,180,1,1,0.005,-0.002,0.00001,0,0,1,0,0,0\n";
      TempDir d; auto h = Load(d, ds);
      ScriptedEngine eng{0.2, 0.25, 0.26};
      G4Random::setTheEngine(&eng);
      const G4ThreeVector centre(10 * CLHEP::mm, 20 * CLHEP::mm, 30 * CLHEP::mm);
      CHECK_VEC(h.SampleExitPosition(1, 0, 0, 180, centre, TH, PH), G4ThreeVector(10, 21, 30) * CLHEP::mm, 1e-9);
      CHECK_VEC(h.SampleExitPosition(1, 0, 0, 180, centre, TH, PH), G4ThreeVector(10, 21, 30) * CLHEP::mm, 1e-9);
      CHECK_VEC(h.SampleExitPosition(1, 0, 0, 180, centre, TH, PH), G4ThreeVector(10, 18, 30.01) * CLHEP::mm, 1e-9);
    }},
    {"err_missing_files", [] {
      Dataset a = hfssfix::Mini500(); a.ff1.clear();
      TempDir d1; ExpectG4Exception("BBR001", [&] { Load(d1, a); }, "BBRHFSSData");
      Dataset b = hfssfix::Mini500(); b.wg0.clear();
      TempDir d2; ExpectG4Exception("BBR002", [&] { Load(d2, b); }, "BBRHFSSData");
    }},
    {"err_no_keys", [] {
      Dataset ds; ds.ff0 = hfssfix::FF_HDR; ds.ff1 = hfssfix::FF_HDR; ds.wg0 = hfssfix::WG_HDR; ds.wg1 = hfssfix::WG_HDR;
      TempDir d; ExpectG4Exception("BBR000", [&] { Load(d, ds); }, "BBRHFSSData");
    }},
    {"err_key_mismatch", [] {
      Dataset ds = hfssfix::Mini500();
      // Both waveguide files use key (0,135), both far fields (0,180). Changing
      // only wg0 would make the Ephi=1 waveguide row unpaired, which is BBR012.
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,135,3,4,0,0,0,0,1,0,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,135,1,4,0,0,0,0,0,1,0,0,0\n";
      TempDir d; ExpectG4Exception("BBR007", [&] { Load(d, ds); }, "BBRHFSSData");
    }},
    {"err_freq_column", [] {
      Dataset ds = hfssfix::Mini500();
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n501GHz,1,0,180,1,4,0,0.001,0,0,0,1,0,0,0\n";
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,3,4,0,0.001,0,0,1,0,0,0,0\n";
      TempDir d; ExpectG4Exception("BBR009", [&] { Load(d, ds); }, "BBRHFSSData");   // row 2 (D9: every row is checked)
      Dataset ok = hfssfix::Mini500();
      ok.wg0 = std::string(hfssfix::WG_HDR) + "500.5GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n";   // 0.1 %: accepted
      TempDir d2; auto h = Load(d2, ok); CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 0.75, 1e-15);
    }},
    {"err_ephi1_row_mismatch", [] {
      // D9: the Ephi=1 exit point at index 0 sits at a different Y than the Ephi=0 one.
      Dataset ds = hfssfix::Mini500();
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0.001,0,0,0,1,0,0,0\n";
      TempDir d; ExpectG4Exception("BBR012", [&] { Load(d, ds); }, "BBRHFSSData");
      // and a far-field grid mismatch
      Dataset fs = hfssfix::Mini500();
      fs.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,30,90,1,0,0,0\n";
      TempDir d2; ExpectG4Exception("BBR012", [&] { Load(d2, fs); }, "BBRHFSSData");
    }},
    {"err_bad_number", [] {
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,abc,4,0,0,0,0,1,0,0,0,0\n";
      TempDir d; ExpectG4Exception("BBR013", [&] { Load(d, ds); }, "BBRHFSSData");
      Dataset nn = hfssfix::Mini500();
      nn.ff0 = std::string(hfssfix::FF_HDR) + "500GHz,0,0,180,0,90,0,0,nan,0\n";
      TempDir d2; ExpectG4Exception("BBR013", [&] { Load(d2, nn); }, "BBRHFSSData");
      // CDX-001: an infinite amplitude or power, of either sign, is refused at load too.
      for (const char* inf : {"inf", "-inf", "INF"}) {
        Dataset fi = hfssfix::Mini500();
        fi.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90," + inf + ",0,0,0\n";
        TempDir d4; ExpectG4Exception("BBR013", [&] { Load(d4, fi); }, "BBRHFSSData", "rEphi_real");
        Dataset wi = hfssfix::Mini500();
        wi.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180," + inf + ",4,0,0,0,0,1,0,0,0,0\n";
        TempDir d5; ExpectG4Exception("BBR013", [&] { Load(d5, wi); }, "BBRHFSSData", "OutgoingPower");
      }
      Dataset sr = hfssfix::Mini500();
      sr.wg0 += "500GHz,0,0,180,3,4,0,0.001,0\n";   // a short row: 9 of 15 fields
      TempDir d3; ExpectG4Exception("BBR013", [&] { Load(d3, sr); }, "BBRHFSSData");
    }},
    {"zero_support_fatal", [] {
      // CDX-008: a transmitting state must have outward far-field intensity and exit field.
      Dataset ff0 = hfssfix::Mini500();   // T0 = 3/4 with an all-zero Ephi=0 far field
      ff0.ff0 = std::string(hfssfix::FF_HDR) + "500GHz,0,0,180,0,90,0,0,0,0\n";
      TempDir d1; ExpectG4Exception("BBR027", [&] { Load(d1, ff0); }, "BBRHFSSData",
                                    "key (0, 180): T0 = 0.75 > 0, but its far field has no outward intensity");
      Dataset plane = hfssfix::Mini500();   // the only far-field row lies in the exit-face plane (Phi = -90)
      plane.ff0 = ReplaceAll(plane.ff0, ",0,180,0,90,", ",0,180,-90,90,");
      plane.ff1 = ReplaceAll(plane.ff1, ",0,180,0,90,", ",0,180,-90,90,");
      TempDir d2; ExpectG4Exception("BBR027", [&] { Load(d2, plane); }, "BBRHFSSData", "far field has no outward intensity");
      Dataset ex1 = hfssfix::Mini500();   // T1 = 1/4 with an all-zero Ephi=1 exit field
      ex1.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,0,0,0,0\n";
      TempDir d3; ExpectG4Exception("BBR027", [&] { Load(d3, ex1); }, "BBRHFSSData",
                                    "key (0, 180): T1 = 0.25 > 0, but its exit field is zero at every exit point");
      // A zero-filled point outside the section is never a fallback: with no field anywhere, the load fails.
      Dataset outside = hfssfix::Mini500();
      outside.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,0,0,0,0,0,0,0,0,0\n500GHz,0,0,180,3,4,0,0.01,0,0,0,0,0,0,0\n";
      outside.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,1,4,0,0.01,0,0,0,0,0,0,0\n";
      TempDir d4; ExpectG4Exception("BBR027", [&] { Load(d4, outside); }, "BBRHFSSData", "T0 = 0.75 > 0, but its exit field is zero");
      // A state that does not transmit needs no support: T1 = 0 with no Ephi=1 field loads.
      Dataset dark = hfssfix::Mini500();
      dark.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90,0,0,0,0\n";
      dark.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0,4,0,0,0,0,0,0,0,0,0\n";
      TempDir d5; CHECK_NEAR(Load(d5, dark).GetTransmittance(0, 1, 0, 180), 0., 0.);
      // Mixed-state cancellation: F1 = -F0 and E1 = -E0 with T0 = 3/4, T1 = 1/4, so at E_theta = E_phi
      // T = (sqrt(3/4) - sqrt(1/4))^2 / 2 = 0.067 > 0 while both samplers have zero weight everywhere.
      Dataset cancel = hfssfix::Mini500();
      cancel.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90,0,0,-1,0\n";
      cancel.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,-1,0,0,0,0\n";
      TempDir d6; const auto h = Load(d6, cancel);
      const double c = std::sqrt(0.5);
      CHECK_NEAR(h.GetTransmittance(c, c, 0, 180), 0.5 * std::pow(std::sqrt(0.75) - 0.5, 2), 1e-12);
      G4ThreeVector pol;
      ExpectG4Exception("BBR027", [&] { h.SampleOutgoingDirection(c, c, 0, 180, G4ThreeVector(0, 0, 1), G4ThreeVector(0, 1, 0),
                                                                  G4ThreeVector(1, 0, 0), pol); },
                        "BBRHFSSData::SampleOutgoingDirection", "the basis far fields cancel");
      ExpectG4Exception("BBR027", [&] { h.SampleExitPosition(c, c, 0, 180, G4ThreeVector(), G4ThreeVector(0, 1, 0),
                                                             G4ThreeVector(0, 0, 1)); },
                        "BBRHFSSData::SampleExitPosition", "the basis exit fields cancel at every point");
    }},
    {"real_data_smoke", [] {
      // Independent numbers: awk over the shipped CSVs (first row per key, Out/In). The
      // load divides (0,180) by its largest coherent T, 1.05451 (raw T0 above 1).
      const char* root = std::getenv("BBRSIM_TEST_DATA");
      if (!root) { CHECK(false && "BBRSIM_TEST_DATA not set"); return; }
      BBRHFSSData h(std::string(root) + "/waveguides", "InfParallelPlate_crack1Rohan_500GHz", 500.);
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 1.0, 1e-9);          // raw 1.05451 normalized (1 - 2e-13)
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 180), 1.29671e-10, 1e-14);  // raw 1.36739e-10 / 1.05451
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 135), 0.759519, 1e-6);
      CHECK_NEAR(h.GetTransmittance(0, 1, 90, 135), 0.773205, 1e-6);
    }},
    {"transmittance_coherent", [] {
      // Identical exit fields for both polarizations (rho = 1): the transmitted wave is the
      // projection of the polarization onto one axis, so T(1,1)/sqrt2 = T0 + T1 and T(1,-1)/sqrt2 = 0.
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,0.4,1,0,0,0,0,1,0,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.4,1,0,0,0,0,1,0,0,0,0\n";
      TempDir d; auto h = Load(d, ds);
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 0.4, 0.);   // lambda_max = 0.8 <= 1: untouched
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 180), 0.4, 0.);
      CHECK_NEAR(h.GetTransmittance(R2, R2, 0, 180), 0.8, 1e-12);
      CHECK_NEAR(h.GetTransmittance(R2, -R2, 0, 180), 0.0, 1e-12);
      // Orthogonal exit fields (Mini500 as written: y vs z): no cross term.
      TempDir d2; auto o = Load(d2, hfssfix::Mini500());
      CHECK_NEAR(o.GetTransmittance(R2, R2, 0, 180), 0.5, 1e-15);
      CHECK_NEAR(o.GetTransmittance(R2, -R2, 0, 180), 0.5, 1e-15);
      // Quadrature phase (E1 = i E0): Re rho = 0, no cross term either.
      Dataset q = hfssfix::Mini500();
      q.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,0.4,1,0,0,0,0,1,0,0,0,0\n";
      q.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.4,1,0,0,0,0,0,0,0,1,0\n";
      TempDir d3; auto qq = Load(d3, q);
      CHECK_NEAR(qq.GetTransmittance(R2, R2, 0, 180), 0.4, 1e-12);
      // Load-time normalization: a key whose largest T over linear polarizations (the top
      // eigenvalue of [[T0, c], [c, T1]], c = sqrt(T0 T1) Re rho) exceeds 1 has T0 and T1
      // divided by it. rho = 1 with T0 = 0.8, T1 = 0.4: the in-phase T would be 1.2 although
      // neither T0 nor T1 exceeds 1; after the division it is 1 and T1/T0 is kept.
      Dataset s = hfssfix::Mini500();
      s.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,0.8,1,0,0,0,0,1,0,0,0,0\n";
      s.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.4,1,0,0,0,0,1,0,0,0,0\n";
      TempDir d4; auto ss = Load(d4, s);
      const double t0 = ss.GetTransmittance(1, 0, 0, 180), t1 = ss.GetTransmittance(0, 1, 0, 180);
      CHECK_NEAR(t0, 0.8 / 1.2, 1e-15);
      CHECK_NEAR(t1 / t0, 0.5, 1e-15);
      CHECK_NEAR(ss.GetTransmittance(std::sqrt(2. / 3.), std::sqrt(1. / 3.), 0, 180), 1.0, 1e-15);   // in phase
      CHECK_NEAR(ss.GetTransmittance(std::sqrt(1. / 3.), -std::sqrt(2. / 3.), 0, 180), 0.0, 1e-15);  // orthogonal
      // rho = 0 with T0 = 1.25, T1 = 0.5: lambda_max = 1.25, so T0 -> 1 and T1 -> 0.4 (the
      // per-polarization cap this replaces left T1 at 0.5).
      Dataset m = hfssfix::Mini500();
      m.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,1.25,1,0,0,0,0,1,0,0,0,0\n";
      m.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,0.5,1,0,0,0,0,0,1,0,0,0\n";
      TempDir d5; auto mm = Load(d5, m);
      CHECK_NEAR(mm.GetTransmittance(1, 0, 0, 180), 1.0, 0.);
      CHECK_NEAR(mm.GetTransmittance(0, 1, 0, 180), 0.4, 0.);
    }},
    {"real_data_diagonal_key", [] {
      // (45,135) of crack1: the two basis exit fields are 99.96 % correlated (survey), so an
      // in-phase diagonal polarization transmits ~T0+T1 and the orthogonal one ~0.
      const char* root = std::getenv("BBRSIM_TEST_DATA");
      if (!root) { CHECK(false && "BBRSIM_TEST_DATA not set"); return; }
      BBRHFSSData h(std::string(root) + "/waveguides", "InfParallelPlate_crack1Rohan_500GHz", 500.);
      const double T0 = h.GetTransmittance(1, 0, 45, 135), T1 = h.GetTransmittance(0, 1, 45, 135);
      CHECK_NEAR(T0, 0.4073, 5e-4); CHECK_NEAR(T1, 0.3939, 5e-4);
      CHECK_NEAR(h.GetTransmittance(R2, R2, 45, 135), T0 + T1, 0.002);      // 0.80
      CHECK_NEAR(h.GetTransmittance(R2, -R2, 45, 135), 0.0, 0.001);         // 2.3e-4 (real data)
      // and unchanged at the orthogonal keys
      CHECK_NEAR(h.GetTransmittance(R2, R2, 0, 180), 0.5, 0.01);
      CHECK_NEAR(h.GetTransmittance(R2, R2, 90, 135), 0.5 * 0.773205, 0.01);
      // The load-time normalization in the real data, both cracks. At (45,180) the raw
      // in-phase T exceeds 1 (T0 + T1 = 1.025 / 1.006, Re rho ~ 1) while neither T0 nor T1
      // does; T0 and T1 are divided by lambda_max (1.02510 / 1.00570), so the random-
      // polarization mean (the psi-average over linear polarizations) is (T0 + T1)/2 with
      // no clamp acting: 0.500094 / 0.500033, where the per-polarization cap gave 0.512645 /
      // 0.502883. At (0,180) the raw T0 alone exceeds 1 (1.05451 / 1.00873): T0 is 1 to 2e-12.
      const std::pair<const char*, double> cracks[] = {
        {"InfParallelPlate_crack1Rohan_500GHz", 0.500094}, {"InfParallelPlate_crack2_500GHz", 0.500033}};
      for (const auto& [stem, mean45] : cracks) {
        BBRHFSSData c(std::string(root) + "/waveguides", stem, 500.);
        const double t0 = c.GetTransmittance(1, 0, 45, 180), t1 = c.GetTransmittance(0, 1, 45, 180);
        double sum = 0.;
        for (int i = 0; i < 360; ++i)
          sum += c.GetTransmittance(std::cos(i * CLHEP::deg), std::sin(i * CLHEP::deg), 45, 180);
        CHECK_NEAR(sum / 360., 0.5 * (t0 + t1), 1e-12);   // the clamp never acts
        CHECK_NEAR(sum / 360., mean45, 1e-6);
        CHECK_NEAR(c.GetTransmittance(1, 0, 0, 180), 1.0, 1e-9);
      }
    }},
    {"sidecar_checks_pass", [] {
      TempDir d1; CHECK_NEAR(LoadChecked(d1, hfssfix::Mini500()).GetFrequencyGHz(), 500., 0);
      TempDir d2; LoadChecked(d2, hfssfix::NineKeyGrid());             // 9 far-field rows, 4 exit points per key
      TempDir d3; LoadChecked(d3, hfssfix::WrapKeys("500", "0.01"));   // the section grows to the 10 mm rows
    }},
    {"sidecar_c1_header", [] {
      Dataset ds = hfssfix::Mini500();
      const auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      ds.ff0.replace(ds.ff0.find(",Theta,"), 7, ",Thet,");
      TempDir d; ExpectG4Exception("BBR013", [&] { LoadChecked(d, ds, &p); }, "BBRHFSSData",
                                   "c_500GHz_Ephi=0/far_field.csv differs from the columns its sidecar test declares: "
                                   "column 5 is 'Thet', the sidecar declares 'Theta'");
      // waveguide.csv goes through C1 as well: a renamed column, then one column too many.
      Dataset wg = hfssfix::Mini500();
      wg.wg1 = ReplaceAll(wg.wg1, ",Ez_imag\n", ",Ez_imga\n");
      TempDir d2; ExpectG4Exception("BBR013", [&] { LoadChecked(d2, wg, &p); }, "BBRHFSSData",
                                    "c_500GHz_Ephi=1/waveguide.csv differs from the columns its sidecar test declares: "
                                    "column 14 is 'Ez_imga', the sidecar declares 'Ez_imag'");
      wg = hfssfix::Mini500();
      wg.wg0 = ReplaceAll(wg.wg0, ",Ez_imag\n", ",Ez_imag,Extra\n");
      TempDir d3; ExpectG4Exception("BBR013", [&] { LoadChecked(d3, wg, &p); }, "BBRHFSSData",
                                    "c_500GHz_Ephi=0/waveguide.csv differs from the columns its sidecar test declares: "
                                    "the header has 16 columns, the sidecar declares 15");
    }},
    {"sidecar_c1_bom", [] {
      // A UTF-8 byte-order mark before the header and CRLF line ends, as a Windows tool writes
      // them: pandas strips the BOM, so the Python check passes such a file and C1 must too.
      Dataset ds = hfssfix::Mini500();
      const auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      for (std::string* csv : {&ds.ff0, &ds.ff1, &ds.wg0, &ds.wg1}) *csv = "\xEF\xBB\xBF" + ReplaceAll(*csv, "\n", "\r\n");
      TempDir d; CHECK_NEAR(LoadChecked(d, ds, &p).GetFrequencyGHz(), 500., 0);
      // Only the mark is dropped: a header that differs after it still fails.
      ds.ff1 = ReplaceAll(ds.ff1, ",Phi,", ",Phj,");
      TempDir d2; ExpectG4Exception("BBR013", [&] { LoadChecked(d2, ds, &p); }, "BBRHFSSData",
                                    "c_500GHz_Ephi=1/far_field.csv differs from the columns its sidecar test declares: "
                                    "column 4 is 'Phj', the sidecar declares 'Phi'");
    }},
    {"sidecar_c2_key_off_grid", [] {
      auto p = hfssfix::SidecarFrom("c_500GHz", hfssfix::Mini500());
      p.thetas = {135.};
      TempDir d; ExpectG4Exception("BBR007", [&] { LoadChecked(d, hfssfix::Mini500(), &p); },
                                   "BBRHFSSData", "incidence key (0, 180) is not on the declared");
    }},
    {"sidecar_c2_key_missing", [] {
      // CDX-021: ThreeKey's sidecar declares {0, 90} x {135, 180}, but (90,135) has no rows; the
      // nearest-key lookup would serve (90,180)'s or (0,135)'s tables for it.
      TempDir d; ExpectG4Exception("BBR007", [&] { LoadChecked(d, hfssfix::ThreeKey()); }, "BBRHFSSData",
                                   "1 declared incidence key(s) have no rows: (90, 135)");
      TempDir d2; LoadChecked(d2, hfssfix::FourKey());   // the same grid, complete
    }},
    {"sidecar_c3_far_field", [] {
      const Dataset ds = hfssfix::NineKey();
      const auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      auto q = p; q.thMax = 100.;   // the Theta = 110 rows lie outside
      TempDir d1; ExpectG4Exception("BBR012", [&] { LoadChecked(d1, ds, &q); },
                                    "BBRHFSSData", "(Phi, Theta) = (60, 110) of key (0, 90) lies outside the declared ranges "
                                                   "Phi in [-40, 60], Theta in [40, 100] deg");
      // A value just past the bound prints in full, not rounded onto the bound.
      Dataset past = ds;
      past.ff0 = ReplaceAll(past.ff0, ",60,110,", ",60,110.000001,");
      past.ff1 = ReplaceAll(past.ff1, ",60,110,", ",60,110.000001,");
      TempDir d4; ExpectG4Exception("BBR012", [&] { LoadChecked(d4, past, &p); }, "BBRHFSSData",
                                    "(Phi, Theta) = (60, 110.000001) of key (0, 90) lies outside the declared ranges "
                                    "Phi in [-40, 60], Theta in [40, 110] deg");
      q = p; q.ffPerKey = 2;
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadChecked(d2, ds, &q); },
                                    "BBRHFSSData", "has 3 far-field rows; points_per_key is 2");
      q = p; q.phCount = 4;
      TempDir d3; ExpectG4Exception("BBR012", [&] { LoadChecked(d3, ds, &q); },
                                    "BBRHFSSData", "the far field has 3 Phi and 3 Theta values; the sidecar declares 4 and 3");
    }},
    {"sidecar_c4_exit_rows", [] {
      const Dataset ds = hfssfix::NineKey();
      auto q = hfssfix::SidecarFrom("c_500GHz", ds);
      q.exitPerKey = 1;
      TempDir d1; ExpectG4Exception("BBR012", [&] { LoadChecked(d1, ds, &q); },
                                    "BBRHFSSData", "has 2 exit points; points_per_key_retained is 1");
      Dataset off = hfssfix::Mini500();   // X = 1 um: not on the exit plane
      off.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,1e-06,0,0,0,1,0,0,0,0\n";
      off.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,1e-06,0,0,0,0,1,0,0,0\n";
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadChecked(d2, off); },
                                    "BBRHFSSData", "has X = 1e-06 m; the sidecar declares x_e = 0");
      auto r = hfssfix::SidecarFrom("c_500GHz", ds);
      r.yMax = 0.0005;              // the Y = 1 mm points lie beyond the declared grid
      TempDir d3; ExpectG4Exception("BBR012", [&] { LoadChecked(d3, ds, &r); },
                                    "BBRHFSSData", "(Y, Z) = (0.001, 1e-05) m of key (0, 90) lies outside the declared grid "
                                                   "Y in [-0.002, 0.0005], Z in [5e-06, 1e-05] m");
      Dataset past = ds;   // 1 nm past the declared Y = 1 mm edge: printed in full
      past.wg0 = ReplaceAll(past.wg0, ",0,0.001,1e-05,", ",0,0.001000001,1e-05,");
      past.wg1 = ReplaceAll(past.wg1, ",0,0.001,1e-05,", ",0,0.001000001,1e-05,");
      const auto declared = hfssfix::SidecarFrom("c_500GHz", ds);   // Y in [-2, 1] mm, as NineKey writes it
      TempDir d5; ExpectG4Exception("BBR012", [&] { LoadChecked(d5, past, &declared); }, "BBRHFSSData",
                                    "(Y, Z) = (0.001000001, 1e-05) m of key (0, 90) lies outside the declared grid "
                                    "Y in [-0.002, 0.001], Z in [5e-06, 1e-05] m");
      r = hfssfix::SidecarFrom("c_500GHz", ds);
      r.zCount = 3;
      TempDir d4; ExpectG4Exception("BBR012", [&] { LoadChecked(d4, ds, &r); },
                                    "BBRHFSSData", "the exit grid has 2 Y and 2 Z values; the sidecar declares 2 and 3");
      // With outside_points "omitted" the producer declares the whole export lattice and drops
      // the lattice points outside the section, so a round face whose radius is not a multiple
      // of the step leaves its outermost lattice columns without a retained point: the retained
      // distinct counts may be smaller than the declared ones, never larger. Here R = 15 um at a
      // 10 um step: lattice Y in [-20, 20] um, 5 columns declared, 3 retained (-10, 0, +10 um).
      Dataset disc = hfssfix::Mini500();
      disc.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,0,-1e-05,0,0,1,0,0,0,0\n"
                 "500GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,3,4,0,1e-05,0,0,1,0,0,0,0\n";
      disc.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,-1e-05,0,0,0,1,0,0,0\n"
                 "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,1,4,0,1e-05,0,0,0,1,0,0,0\n";
      auto w = hfssfix::SidecarFrom("c_500GHz", disc);   // Y in [-10, 10] um, 3 columns, 3 points per key
      w.disc = true; w.radius = 1.5e-5;
      w.yMin = -2e-5; w.yMax = 2e-5; w.yCount = 5;
      w.outsidePoints = "omitted";
      TempDir d6; LoadChecked(d6, disc, &w);              // accepted: 3 retained of 5 declared
      w.outsidePoints = "none";                            // every lattice point written: equality
      TempDir d7; ExpectG4Exception("BBR012", [&] { LoadChecked(d7, disc, &w); },
                                    "BBRHFSSData", "the exit grid has 3 Y and 1 Z values; the sidecar declares 5 and 1");
      w.outsidePoints = "zero";
      TempDir d8; ExpectG4Exception("BBR012", [&] { LoadChecked(d8, disc, &w); },
                                    "BBRHFSSData", "the exit grid has 3 Y and 1 Z values; the sidecar declares 5 and 1");
      // More retained values than the declared lattice has columns is wrong under "omitted" too.
      w.outsidePoints = "omitted"; w.yCount = 2;           // a lattice of two columns at +-20 um
      TempDir d9; ExpectG4Exception("BBR012", [&] { LoadChecked(d9, disc, &w); },
                                    "BBRHFSSData", "the exit grid has 3 Y and 1 Z values; the sidecar declares 2 and 1 "
                                                   "(with outside_points omitted the retained counts may be smaller, not larger)");
    }},
    {"sidecar_c3_c4_per_key", [] {
      // CDX-022: the dataset-wide counts pass each of these; the per-key checks do not.
      const std::vector<std::string> one = {"0.001,0"};
      TempDir d1; ExpectG4Exception("BBR012", [&] { LoadChecked(d1, Rows({"0,60", "0,60", "30,60", "30,120"}, one)); },
                                    "BBRHFSSData", "far-field point (Phi, Theta) = (0, 60) appears twice in key (0, 180)");
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadChecked(d2, Rows({"0,60", "0,120", "0,61", "30,60", "30,120", "30,61"}, one)); },
                                    "BBRHFSSData", "far-field row (Phi, Theta) = (0, 61) of key (0, 180) is not a node of the declared "
                                                   "grid Phi 0 to 30 in 2 nodes x Theta 60 to 120 in 3 nodes");
      // Two keys whose far fields differ in compensating ways: the union is the 2 x 2 grid.
      TempDir d3; ExpectG4Exception("BBR012", [&] {
        LoadChecked(d3, Keys({{"0,135", {"0,60", "30,120"}}, {"0,180", {"30,60", "0,120"}}}, one)); },
        "BBRHFSSData", "key (0, 135) holds 2 of the 4 nodes of the declared far-field grid");
      const std::vector<std::string> ff = {"0,90"};
      TempDir d4; ExpectG4Exception("BBR012", [&] {
        LoadChecked(d4, Rows(ff, {"-0.001,-5e-06", "-0.001,5e-06", "0,-5e-06", "0,-5e-06", "0.001,-5e-06", "0.001,5e-06"})); },
        "BBRHFSSData", "exit point (Y, Z) = (0, -5e-06) m appears twice in key (0, 180)");
      TempDir d5; ExpectG4Exception("BBR012", [&] { LoadChecked(d5, Rows(ff, {"-0.001,0", "0,0", "0.0015,0"})); },
                                    "BBRHFSSData", "exit point (Y, Z) = (0, 0) m of key (0, 180) is not a node of the declared lattice");
      // Valid: the pole rows Theta = 0 at two Phi (one direction, two grid nodes), and single-node axes.
      TempDir d6; LoadChecked(d6, Rows({"0,0", "90,0", "0,90", "90,90"}, one));
      // "omitted": a lattice node inside the section without a retained point is missing.
      Dataset disc = Rows(ff, {"-1e-05,0", "1e-05,0"});
      auto w = hfssfix::SidecarFrom("c_500GHz", disc);
      w.disc = true; w.radius = 1.5e-5; w.yMin = -2e-5; w.yMax = 2e-5; w.yCount = 5; w.outsidePoints = "omitted";
      TempDir d7; ExpectG4Exception("BBR012", [&] { LoadChecked(d7, disc, &w); }, "BBRHFSSData",
                                    "key (0, 180) holds 2 exit points; every key holds each of the 3 declared lattice nodes "
                                    "inside the section once");
    }},
    {"sidecar_c5_outside_section", [] {
      const Dataset ds = hfssfix::WrapKeys("500", "0.01");   // exit points at Y = 10 mm
      auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      p.yHalf = 5e-3;                                       // the declared section ends at 5 mm
      TempDir d; ExpectG4Exception("BBR025", [&] { LoadChecked(d, ds, &p); }, "BBRHFSSData",
                                   "6 exit point(s) lie outside the declared cross-section, e.g. (Y, Z) = (0.01, 1e-05) m");
      // Two offenders in one key: the message names the first, as the Python check does.
      Dataset two = hfssfix::Mini500();
      two.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n"
                "500GHz,0,0,180,3,4,0,0.01,0,0,1,0,0,0,0\n500GHz,0,0,180,3,4,0,0.02,0,0,1,0,0,0,0\n";
      two.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n"
                "500GHz,1,0,180,1,4,0,0.01,0,0,0,1,0,0,0\n500GHz,1,0,180,1,4,0,0.02,0,0,0,1,0,0,0\n";
      auto q = hfssfix::SidecarFrom("c_500GHz", two);
      q.yHalf = 4.5e-3;
      TempDir d2; ExpectG4Exception("BBR025", [&] { LoadChecked(d2, two, &q); }, "BBRHFSSData",
                                    "2 exit point(s) lie outside the declared cross-section, e.g. (Y, Z) = (0.01, 0) m");
    }},
    {"sidecar_c5_zero_outside", [] {
      // An exit point outside the section with no field at all: exempt only when outside_points is "zero".
      Dataset ds = hfssfix::Mini500();
      ds.wg0 = std::string(hfssfix::WG_HDR) + "500GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,3,4,0,0.01,0,0,0,0,0,0,0\n";
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,1,4,0,0.01,0,0,0,0,0,0,0\n";
      auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      p.yHalf = 4.5e-3;                                    // the point at Y = 10 mm lies outside
      p.outsidePoints = "zero";
      TempDir d1; LoadChecked(d1, ds, &p);                 // exempt: all six components zero in both polarizations
      p.outsidePoints = "omitted";
      TempDir d2; ExpectG4Exception("BBR025", [&] { LoadChecked(d2, ds, &p); },
                                    "BBRHFSSData", "1 exit point(s) lie outside the declared cross-section, e.g. (Y, Z) = (0.01, 0) m");
      // "zero" exempts only a point with no field in EITHER polarization: one Ephi=1 component is enough to refuse it.
      Dataset one = ds;
      one.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,1,4,0,0.01,0,0,0,0,0,0,1e-30\n";
      p.outsidePoints = "zero";
      TempDir d3; ExpectG4Exception("BBR025", [&] { LoadChecked(d3, one, &p); },
                                    "BBRHFSSData", "1 exit point(s) lie outside");
    }},
    {"exit_clamp_rectangle", [] {
      // OR-020 (C2, C3): an exit point on a wall of the declared section, or within C5's tolerance
      // outside it, is moved to 1 nm inside, per axis; every other point is exactly unchanged.
      const double yw = 4.5e-3, zw = 2.5e-5, d = 1e-9;
      CHECK(BBRCrossSection::kClampInset_m == d);           // the 1 nm of the log line (OR-020)
      std::vector<YZ> pts;                                   // Y outer, Z inner: row = 3 iy + iz
      for (const double y : {-4.5e-3, -2.25e-3, 0., 2.25e-3, 4.5e-3})
        for (const double z : {-2.5e-5, 0., 2.5e-5}) pts.push_back({y, z});
      const std::size_t M = pts.size();
      CHECK(M == 15);
      const Dataset ds = UniformLattice(pts);
      const auto clamp1 = [&](double v, double half) { return std::abs(v) > half - d ? std::copysign(half - d, v) : v; };
      ScriptedEngine eng;
      G4Random::setTheEngine(&eng);

      auto p = hfssfix::SidecarFrom("c_500GHz", ds);        // widens the section to the lattice: set it after
      p.yHalf = yw; p.zHalf = zw;
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir dir; auto h = LoadChecked(dir, ds, &p);
      G4iosSetDestination(nullptr);
      CHECK(Count(cap.out, "clamped") == 1);
      CHECK(cap.out.find(ClampLine("c_500GHz", 12)) != std::string::npos);

      std::size_t moved = 0;
      for (std::size_t j = 0; j < M; ++j) {
        const auto r = DrawRow(h, eng, j, M);
        const double y = pts[j].first, z = pts[j].second;
        if (std::abs(y) == yw || std::abs(z) == zw) {
          ++moved;
          CHECK_NEAR(r.y(), clamp1(y, yw) * CLHEP::m, 1e-11);
          CHECK_NEAR(r.z(), clamp1(z, zw) * CLHEP::m, 1e-11);
        } else {
          CHECK(r.y() == y * CLHEP::m);
          CHECK(r.z() == z * CLHEP::m);
        }
      }
      CHECK(moved == 12);
      // The corner moves on both axes, a long-wall point on one, the three points between the walls not at all.
      const auto corner = DrawRow(h, eng, 14, M);           // (+4.5 mm, +25 um)
      CHECK_NEAR(corner.y(), (yw - d) * CLHEP::m, 1e-11);
      CHECK_NEAR(corner.z(), (zw - d) * CLHEP::m, 1e-11);
      const auto end = DrawRow(h, eng, 13, M);              // (+4.5 mm, 0)
      CHECK_NEAR(end.y(), (yw - d) * CLHEP::m, 1e-11);
      CHECK(end.z() == 0.);
      for (const std::size_t j : {std::size_t(4), std::size_t(7), std::size_t(10)}) {   // (-2.25 mm, 0), (0, 0), (+2.25 mm, 0)
        const auto r = DrawRow(h, eng, j, M);
        CHECK(r.y() == pts[j].first * CLHEP::m);
        CHECK(r.z() == 0.);
      }

      // Past the wall, inside C5's relative 1e-6: the end columns at 4.5 mm lie 5e-7 beyond a section of
      // y_half = 4.5 mm / (1 + 5e-7) and still land 1 nm inside it.
      const double yw2 = 4.5e-3 / (1. + 5e-7);
      p.yHalf = yw2;
      Capture cap2;
      G4iosSetDestination(&cap2);
      TempDir dir2; auto h2 = LoadChecked(dir2, ds, &p);
      G4iosSetDestination(nullptr);
      CHECK(cap2.out.find(ClampLine("c_500GHz", 12)) != std::string::npos);
      for (std::size_t j = 0; j < M; ++j) {
        const auto r = DrawRow(h2, eng, j, M);
        CHECK_NEAR(r.y(), clamp1(pts[j].first, yw2) * CLHEP::m, 1e-11);
        CHECK_NEAR(r.z(), clamp1(pts[j].second, zw) * CLHEP::m, 1e-11);
      }
      const auto end2 = DrawRow(h2, eng, 13, M);            // (+4.5 mm, 0)
      CHECK_NEAR(end2.y(), (yw2 - d) * CLHEP::m, 1e-11);
      CHECK(end2.z() == 0.);

      // The band (half - 1 nm, half): a point 0.5 nm inside a wall moves out to half - 1 nm, a point exactly
      // 1 nm inside stays, and so do the columns between. Complete lattice, 5 Y columns by 3 Z columns, from
      // 1 nm inside the low wall to 0.5 nm inside the high one, on both axes; the sidecar declares yw and zw.
      const double yLo = -(yw - d), yHi = yw - 0.5e-9, zLo = -(zw - d), zHi = zw - 0.5e-9;
      const auto node = [](double lo, double hi, int i, int n) { return i == n - 1 ? hi : lo + i * (hi - lo) / (n - 1); };
      std::vector<YZ> band;
      for (int iy = 0; iy < 5; ++iy)
        for (int iz = 0; iz < 3; ++iz) band.push_back({node(yLo, yHi, iy, 5), node(zLo, zHi, iz, 3)});
      const std::size_t Mb = band.size();
      const Dataset dsb = UniformLattice(band);
      auto pb = hfssfix::SidecarFrom("c_500GHz", dsb);
      pb.yHalf = yw; pb.zHalf = zw;
      Capture cap3;
      G4iosSetDestination(&cap3);
      TempDir dir3; auto h3 = LoadChecked(dir3, dsb, &pb);
      G4iosSetDestination(nullptr);
      CHECK(cap3.out.find(ClampLine("c_500GHz", 7)) != std::string::npos);   // Y = high column (3) + Z = high column (5) - the corner
      for (std::size_t j = 0; j < Mb; ++j) {
        const auto r = DrawRow(h3, eng, j, Mb);
        const int iy = int(j) / 3, iz = int(j) % 3;
        if (iy == 4) CHECK_NEAR(r.y(), (yw - d) * CLHEP::m, 1e-11);   // 0.5 nm inside: out to 1 nm inside
        else         CHECK(r.y() == band[j].first * CLHEP::m);        // iy = 0 is exactly 1 nm inside: unchanged
        if (iz == 2) CHECK_NEAR(r.z(), (zw - d) * CLHEP::m, 1e-11);
        else         CHECK(r.z() == band[j].second * CLHEP::m);
      }
      CHECK(eng.Remaining() == 0);
    }},
    {"exit_clamp_disc", [] {
      // The four rim points of a disc land 1 nm inside it, on their own axis; the centre stays.
      const double R = 5e-5, d = 1e-9;
      const std::vector<YZ> pts = {{-5e-5, 0.}, {0., -5e-5}, {0., 0.}, {0., 5e-5}, {5e-5, 0.}};   // the lattice {-R, 0, R}^2 inside the disc
      const std::size_t M = pts.size();
      const Dataset ds = UniformLattice(pts);
      ScriptedEngine eng;
      G4Random::setTheEngine(&eng);

      auto p = hfssfix::SidecarFrom("c_500GHz", ds);        // lattice {-R, 0, R} on both axes, 5 retained points
      p.disc = true; p.radius = R; p.outsidePoints = "omitted";
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir dir; auto h = LoadChecked(dir, ds, &p);
      G4iosSetDestination(nullptr);
      CHECK(Count(cap.out, "clamped") == 1);
      CHECK(cap.out.find(ClampLine("c_500GHz", 4)) != std::string::npos);
      for (std::size_t j = 0; j < M; ++j) {
        const auto r = DrawRow(h, eng, j, M);
        if (j == 2) {
          CHECK(r.y() == 0.);
          CHECK(r.z() == 0.);
        } else {
          CHECK_NEAR(r.y(), std::copysign(pts[j].first != 0. ? R - d : 0., pts[j].first) * CLHEP::m, 1e-11);
          CHECK_NEAR(r.z(), std::copysign(pts[j].second != 0. ? R - d : 0., pts[j].second) * CLHEP::m, 1e-11);
          CHECK_NEAR(std::hypot(r.y(), r.z()), (R - d) * CLHEP::m, 1e-11);
        }
      }

      // A radius that the lattice's rim exceeds by 4e-7 (inside C5's relative 1e-6): the rim is clamped to R2 - 1 nm.
      const double R2 = R / (1. + 4e-7);
      p.radius = R2;
      Capture cap2;
      G4iosSetDestination(&cap2);
      TempDir dir2; auto h2 = LoadChecked(dir2, ds, &p);
      G4iosSetDestination(nullptr);
      CHECK(cap2.out.find(ClampLine("c_500GHz", 4)) != std::string::npos);
      const auto rim = DrawRow(h2, eng, 4, M);              // (+R, 0)
      CHECK_NEAR(rim.y(), (R2 - d) * CLHEP::m, 1e-11);
      CHECK_NEAR(rim.z(), 0., 1e-11);
      const auto mid = DrawRow(h2, eng, 2, M);              // (0, 0)
      CHECK(mid.y() == 0.);
      CHECK(mid.z() == 0.);
      CHECK(eng.Remaining() == 0);
    }},
    {"exit_clamp_disc_off_axis", [] {
      // The clamp of a disc runs along the radius: a rim point off the axes, (30, 40) um, moves to
      // r = R - 1 nm on its own ray, where a per-axis clamp would leave it (each coordinate is far from R - 1 nm).
      // Lattice: 10 um steps from -50 to 50 um in Y and Z, the 81 nodes inside R = 50 um, 12 of them on the rim.
      const double R = 5e-5, d = 1e-9, step = 1e-5;
      std::vector<YZ> pts;
      std::vector<std::pair<int, int>> ij;
      std::size_t onRim = 0;
      for (int iy = -5; iy <= 5; ++iy)
        for (int iz = -5; iz <= 5; ++iz)
          if (iy * iy + iz * iz <= 25) {
            pts.push_back({iy * step, iz * step});
            ij.push_back({iy, iz});
            if (iy * iy + iz * iz == 25) ++onRim;
          }
      const std::size_t M = pts.size();
      CHECK(M == 81);
      CHECK(onRim == 12);
      const auto row = [&](int iy, int iz) {
        for (std::size_t j = 0; j < M; ++j) if (ij[j] == std::make_pair(iy, iz)) return j;
        throw std::runtime_error("no such lattice node");
      };
      const Dataset ds = UniformLattice(pts);
      ScriptedEngine eng;
      G4Random::setTheEngine(&eng);

      auto p = hfssfix::SidecarFrom("c_500GHz", ds);        // 11 x 11 lattice, 81 retained points
      p.disc = true; p.radius = R; p.outsidePoints = "omitted";
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir dir; auto h = LoadChecked(dir, ds, &p);
      G4iosSetDestination(nullptr);
      CHECK(Count(cap.out, "clamped") == 1);
      CHECK(cap.out.find(ClampLine("c_500GHz", 12)) != std::string::npos);

      const double tol = 1e-15 * CLHEP::m;
      std::vector<G4ThreeVector> got;
      for (std::size_t j = 0; j < M; ++j) got.push_back(DrawRow(h, eng, j, M));
      for (std::size_t j = 0; j < M; ++j) {
        const double y = pts[j].first, z = pts[j].second;
        if (ij[j].first * ij[j].first + ij[j].second * ij[j].second == 25) {   // on the rim: along its radius to R - 1 nm
          CHECK_NEAR(got[j].y(), y * (R - d) / R * CLHEP::m, tol);
          CHECK_NEAR(got[j].z(), z * (R - d) / R * CLHEP::m, tol);
          CHECK_NEAR(std::hypot(got[j].y(), got[j].z()), (R - d) * CLHEP::m, tol);
        } else {                                                                // inside: exactly unchanged
          CHECK(got[j].y() == y * CLHEP::m);
          CHECK(got[j].z() == z * CLHEP::m);
        }
      }
      const auto& off = got[row(3, 4)];                      // (30, 40) um
      CHECK_NEAR(off.y(), 3e-5 * (R - d) / R * CLHEP::m, tol);
      CHECK_NEAR(off.z(), 4e-5 * (R - d) / R * CLHEP::m, tol);
      CHECK_NEAR(std::hypot(off.y(), off.z()), (R - d) * CLHEP::m, tol);
      CHECK(off.y() < 3e-5 * CLHEP::m);                      // it moved, and not onto an axis value
      const auto& mid = got[row(1, 2)];                      // (10, 20) um
      CHECK(mid.y() == pts[row(1, 2)].first * CLHEP::m);
      CHECK(mid.z() == pts[row(1, 2)].second * CLHEP::m);
      CHECK(eng.Remaining() == 0);
    }},
    {"exit_clamp_needs_sidecar", [] {
      // The clamp needs the declared section; without a sidecar nothing is moved and nothing is logged.
      std::vector<YZ> pts;
      for (const double y : {-4.5e-3, -2.25e-3, 0., 2.25e-3, 4.5e-3})
        for (const double z : {-2.5e-5, 0., 2.5e-5}) pts.push_back({y, z});
      const std::size_t M = pts.size();
      ScriptedEngine eng;
      G4Random::setTheEngine(&eng);
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir d; auto h = Load(d, UniformLattice(pts));      // direct construction, no sidecar
      G4iosSetDestination(nullptr);
      CHECK(cap.out.find("clamped") == std::string::npos);
      const auto r = DrawRow(h, eng, 13, M);                 // (+4.5 mm, 0)
      CHECK(r.y() == 4.5e-3 * CLHEP::m);
      CHECK(r.z() == 0.);
      const auto c = DrawRow(h, eng, 14, M);                 // (+4.5 mm, +25 um)
      CHECK(c.y() == 4.5e-3 * CLHEP::m);
      CHECK(c.z() == 2.5e-5 * CLHEP::m);
    }},
    {"sidecar_real_data", [] {
      const char* root = std::getenv("BBRSIM_TEST_DATA");
      if (!root) { CHECK(false && "BBRSIM_TEST_DATA not set"); return; }
      const std::string wg = std::string(root) + "/waveguides";
      for (const std::string id : {"InfParallelPlate_crack1Rohan", "InfParallelPlate_crack2"}) {
        const auto sc = BBRDatasetSidecar::Load(wg, id, id + "_500GHz", 500.);
        BBRHFSSData h(wg, id + "_500GHz", 500., &sc);
        CHECK_NEAR(h.GetFrequencyGHz(), 500., 0);
      }
    }},

    // --- schema 2 (OR-034, OR-038, HF-040, HF-041, HF-042) -------------------------------
    {"s2_positions", [] {
      // The 2.0 positions (Ex = 9,10; Ey = 11,12) times the trapezoid weight; the 1.0 map would
      // read Ex = (9, 12), Ey = (10, 13). (0, 180) with E_theta = 1 (its Ephi=1 basis is dark).
      TempDir d; const auto h = LoadV2(d, hfssfix::NineKeyGrid2(), hfssfix::NineKeyParams2(kS2));
      const long n = 200000;
      const P9 p2 = ExitPrediction(1., 0., true);
      P9 p1{};
      for (int i = 0; i < 9; ++i)
        p1[i] = W2(i) * (std::pow(F2(0, i, 0), 2) + std::pow(F2(0, i, 3), 2) + std::pow(F2(0, i, 1), 2) + std::pow(F2(0, i, 4), 2));
      p1 = Normalized(p1);
      const auto c = DrawExits(h, 1., 0., 0., 180., n, 20261019);
      CheckCounts("s2_positions", c, n, p2);
      const double s0 = std::sqrt(p2[0] * (1. - p2[0]) / n);
      std::printf("s2_positions point 0: 2.0 map %.6f, 1.0 map %.6f, separation %.1f sigma\n", p2[0], p1[0],
                  std::abs(p1[0] - p2[0]) / s0);
      CHECK(std::abs(p1[0] - p2[0]) > 10. * s0);                 // the fixture tells the two maps apart
      CHECK(std::abs(double(c[0]) / n - p1[0]) > 10. * s0);      // and the draws reject the 1.0 map
    }},
    {"s2_transmittance", [] {
      // T = OutgoingPower / (IngoingPower * InterceptFactor), kept above 1 and capped per photon.
      TempDir d; const auto h = LoadV2(d, hfssfix::NineKeyGrid2(), hfssfix::NineKeyParams2(kS2));
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 135), std::min(1., 0.6 / std::cos(kPi / 4)), 1e-12);
      CHECK_NEAR(h.GetTransmittance(0, 1, 0, 135), 0.3 / std::cos(kPi / 4), 1e-12);   // no rescale at (0, 135)
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 180), 0.6, 1e-12);
      // (0, 90): c_k = cell_90 (0.0838...), so both linear T exceed 1; the cap holds every mixture at or below 1.
      const double c90 = hfssfix::detail::Cell90F({90., 135., 180.});
      CHECK(0.6 / c90 > 1. && 0.3 / c90 > 1.);
      CHECK(h.GetTransmittance(1, 0, 0, 90) == 1.);
      CHECK(h.GetTransmittance(0, 1, 0, 90) == 1.);   // a load-time rescale would leave T1 below 1
      for (int psi = 0; psi < 360; ++psi) {
        const double T = h.GetTransmittance(std::cos(psi * kPi / 180.), std::sin(psi * kPi / 180.), 0, 90);
        CHECK(T >= 0. && T <= 1.);
      }
    }},
    {"s2_cap_no_rescale", [] {
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir d; const auto h = LoadV2(d, hfssfix::NineKeyGrid2(), hfssfix::NineKeyParams2(kS2));
      G4iosSetDestination(nullptr);
      CHECK(cap.out.find("normalized to 1") == std::string::npos);
      CHECK(Count(cap.out, "clamped") == 1);
    }},
    {"s2_intercept_column", [] {
      const auto p = hfssfix::NineKeyParams2(kS2);
      const auto scale = [](double f) {
        return [f](std::vector<std::string>& v, int) { v[15] = R17(std::stod(v[15]) * f); };
      };
      // One row (the fifth) of (45, 135) Ephi=0 off by 2e-6: not constant within the key.
      Dataset one = hfssfix::NineKeyGrid2();
      one.wg0 = EditKey(one.wg0, 45, 135, [](std::vector<std::string>& v, int r) {
        if (r == 4) v[15] = R17(std::stod(v[15]) * (1. + 2e-6));
      });
      TempDir d1; ExpectG4Exception("BBR012", [&] { LoadV2(d1, one, p); }, "BBRHFSSData", "key (45, 135) polarization 0");
      TempDir d1b; ExpectG4Exception("BBR012", [&] { LoadV2(d1b, one, p); }, "BBRHFSSData", "not constant within the key");
      // Every row of the key off by 2e-6: constant, but not the sidecar's c_key.
      Dataset all = hfssfix::NineKeyGrid2();
      all.wg1 = EditKey(all.wg1, 45, 135, scale(1. + 2e-6));
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadV2(d2, all, p); }, "BBRHFSSData", "key (45, 135) polarization 1");
      TempDir d2b; ExpectG4Exception("BBR012", [&] { LoadV2(d2b, all, p); }, "BBRHFSSData", "transmittance.c_key");
      // Within 1e-6 (5e-7, on one row and on every row of both polarizations): loads.
      Dataset ok = hfssfix::NineKeyGrid2();
      ok.wg0 = EditKey(ok.wg0, 45, 135, [](std::vector<std::string>& v, int r) {
        if (r == 4) v[15] = R17(std::stod(v[15]) * (1. + 5e-7));
      });
      ok.wg0 = EditKey(ok.wg0, 90, 90, scale(1. + 5e-7));
      ok.wg1 = EditKey(ok.wg1, 90, 90, scale(1. + 5e-7));
      TempDir d3; const auto h = LoadV2(d3, ok, p);
      CHECK_NEAR(h.GetTransmittance(1, 0, 0, 135), 0.6 / std::cos(kPi / 4), 1e-12);
    }},
    {"s2_two_t_sources", [] {
      // t_flux 0.6 in the sidecar, OutgoingPower 0.61e-9 over IngoingPower 1e-9 in the CSV.
      Dataset ds = hfssfix::NineKeyGrid2();
      ds.wg0 = EditKey(ds.wg0, 90, 180, [](std::vector<std::string>& v, int) { v[4] = "6.1e-10"; });
      TempDir d; ExpectG4Exception("BBR012", [&] { LoadV2(d, ds, hfssfix::NineKeyParams2(kS2)); }, "BBRHFSSData",
                                   "key (90, 180) polarization 0");
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadV2(d2, ds, hfssfix::NineKeyParams2(kS2)); }, "BBRHFSSData",
                                    "t_flux / c_k");
    }},
    {"s2_dark_outgoing", [] {
      // OR-038: a dark basis needs OutgoingPower exactly 0 on every row; 1e-25 on the third row is refused.
      Dataset ds = hfssfix::NineKeyGrid2();
      ds.wg1 = EditKey(ds.wg1, 0, 180, [](std::vector<std::string>& v, int r) { if (r == 2) v[4] = "1e-25"; });
      TempDir d; ExpectG4Exception("BBR012", [&] { LoadV2(d, ds, hfssfix::NineKeyParams2(kS2)); }, "BBRHFSSData",
                                   "key (0, 180) polarization 1");
      TempDir d2; ExpectG4Exception("BBR012", [&] { LoadV2(d2, ds, hfssfix::NineKeyParams2(kS2)); }, "BBRHFSSData",
                                    "OutgoingPower = 1e-25");
    }},
    {"s2_dark_masking", [] {
      // (0, 180) Ephi=1 is dark with nonzero stored fields: its mixing coefficient is 0 in both
      // samplers, so a tiny bright share samples the Ephi=0 fields alone (and raises no BBR027:
      // any G4Exception fails this case).
      TempDir d; const auto h = LoadV2(d, hfssfix::NineKeyGrid2(), hfssfix::NineKeyParams2(kS2));
      CHECK(h.GetTransmittance(0, 1, 0, 180) == 0.);
      const double a = 0.01, b = 0.99995;
      CHECK_NEAR(h.GetTransmittance(a, b, 0, 180), a * a * 0.6, 1e-15);
      const long n = 200000;
      const P9 pe = ExitPrediction(1., 0., true), pu = ExitPrediction(a, b, true);
      double sep = 0.;
      for (int i = 0; i < 9; ++i) sep = std::max(sep, std::abs(pu[i] - pe[i]) / std::sqrt(pe[i] * (1. - pe[i]) / n));
      std::printf("s2_dark_masking exit: unmasked prediction %.1f sigma away\n", sep);
      CHECK(sep > 10.);
      CheckCounts("s2_dark_masking exit", DrawExits(h, a, b, 0., 180., n, 4242), n, pe);
      const P9 de = DirectionPrediction(1., 0.), du = DirectionPrediction(a, b);
      sep = 0.;
      for (int r = 0; r < 9; ++r) sep = std::max(sep, std::abs(du[r] - de[r]) / std::sqrt(de[r] * (1. - de[r]) / n));
      std::printf("s2_dark_masking direction: unmasked prediction %.1f sigma away\n", sep);
      CHECK(sep > 10.);
      CheckCounts("s2_dark_masking direction", DrawDirections(h, a, b, 0., 180., n, 4243), n, de);
    }},
    {"s2_rho_weighted", [] {
      // (45, 135) with t_flux 0.3 / 0.15 in the sidecar and the CSV powers, fields unchanged: the
      // mixed T stays below the cap, so it shows rho = sum w E_t0.E_t1* / sqrt(sum w|E_t0|^2 sum w|E_t1|^2).
      Dataset ds = hfssfix::NineKeyGrid2();
      ds.wg0 = EditKey(ds.wg0, 45, 135, [](std::vector<std::string>& v, int) { v[4] = "3e-10"; });
      ds.wg1 = EditKey(ds.wg1, 45, 135, [](std::vector<std::string>& v, int) { v[4] = "1.5e-10"; });
      auto p = hfssfix::NineKeyParams2(kS2);
      p.tflux[4] = {0.3, 0.15, 0., 0.};   // phi-major: (45, 135) is the fifth key
      TempDir d; const auto h = LoadV2(d, ds, p);
      const double c = std::cos(kPi / 4), T0 = 0.3 / c, T1 = 0.15 / c;
      const auto rho = [](bool weighted) {
        std::complex<double> x = 0.;
        double p0 = 0., p1 = 0.;
        for (int i = 0; i < 9; ++i) {
          const double w = weighted ? W2(i) : 1.;
          x += w * (Ex2(0, i) * std::conj(Ex2(1, i)) + Ey2(0, i) * std::conj(Ey2(1, i)));
          p0 += w * (std::norm(Ex2(0, i)) + std::norm(Ey2(0, i)));
          p1 += w * (std::norm(Ex2(1, i)) + std::norm(Ey2(1, i)));
        }
        return x.real() / std::sqrt(p0 * p1);
      };
      const double want = 0.5 * T0 + 0.5 * T1 + std::sqrt(T0 * T1) * rho(true);
      const double unweighted = 0.5 * T0 + 0.5 * T1 + std::sqrt(T0 * T1) * rho(false);
      std::printf("s2_rho_weighted: rho %.12f (unweighted %.12f), T %.10f (unweighted %.10f)\n", rho(true), rho(false),
                  want, unweighted);
      CHECK(want > 0. && want < 1.);
      CHECK(std::abs(want - unweighted) > 1e-12);
      CHECK_NEAR(h.GetTransmittance(1. / std::sqrt(2.), 1. / std::sqrt(2.), 45, 135), want, 1e-12);
    }},
    {"s2_weights_before_clamp", [] {
      // The z = +-25 um rows lie on the declared walls: the load clamps them 1 nm inside, after the
      // trapezoid weights were taken from the unclamped lattice. (+1 mm, 0) keeps 1/2, the corner
      // (+1 mm, +25 um) 1/4. (90, 135), E_phi = 1.
      Capture cap;
      G4iosSetDestination(&cap);
      TempDir d; const auto h = LoadV2(d, hfssfix::NineKeyGrid2(), hfssfix::NineKeyParams2(kS2));
      G4iosSetDestination(nullptr);
      CHECK(cap.out.find(ClampLine(kS2, 54)) != std::string::npos);   // 6 points per key, 9 keys
      CHECK(W2(7) == 0.5 && W2(8) == 0.25);
      const long n = 200000;
      const P9 p = ExitPrediction(0., 1., true), pu = ExitPrediction(0., 1., false);
      const double s8 = std::sqrt(p[8] * (1. - p[8]) / n);
      std::printf("s2_weights_before_clamp corner: weighted %.6f, unweighted %.6f, separation %.1f sigma\n", p[8], pu[8],
                  std::abs(pu[8] - p[8]) / s8);
      CHECK(std::abs(pu[8] - p[8]) > 10. * s8);
      double zmax = 0.;
      CheckCounts("s2_weights_before_clamp", DrawExits(h, 0., 1., 90., 135., n, 777, &zmax), n, p);
      CHECK_NEAR(zmax, (2.5e-5 - 1e-9) * CLHEP::m / CLHEP::mm, 1e-12);   // the wall rows were drawn, clamped
    }},
    {"s2_load_bbr027_density", [] {
      // A lit key whose tangential field (Ex, Ey) is zero everywhere but Ez is not: the 2.x exit
      // density w |E_t|^2 is zero, so the load refuses it (BBR027).
      Dataset ds = hfssfix::NineKeyGrid2();
      ds.wg0 = EditKey(ds.wg0, 45, 135, [](std::vector<std::string>& v, int) { v[9] = v[10] = v[11] = v[12] = "0"; });
      TempDir d; ExpectG4Exception("BBR027", [&] { LoadV2(d, ds, hfssfix::NineKeyParams2(kS2)); }, "BBRHFSSData",
                                   "key (45, 135): T0");
    }},
  });
}
