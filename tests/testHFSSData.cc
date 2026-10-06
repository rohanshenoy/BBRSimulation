#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "BBRHFSSData.hh"
#include "BBRDatasetSidecar.hh"
#include "G4SystemOfUnits.hh"
#include "G4coutDestination.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <cmath>
#include <string>

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

// text with every `from` replaced by `to`; throws if `from` is absent, so a stale edit cannot pass.
std::string ReplaceAll(std::string text, const std::string& from, const std::string& to) {
  if (text.find(from) == std::string::npos) throw std::runtime_error("ReplaceAll: '" + from + "' not in the text");
  for (auto p = text.find(from); p != std::string::npos; p = text.find(from, p + to.size())) text.replace(p, from.size(), to);
  return text;
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
      ds.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,0,90,0,0,0,0\n500GHz,1,0,180,30,90,0,0,0,0\n500GHz,1,0,180,60,90,0,0,0,0\n";
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
      ds.ff1 = std::string(hfssfix::FF_HDR) + "500GHz,1,0,180,-90,90,0,0,0,0\n500GHz,1,0,180,0,0,0,0,0,0\n500GHz,1,0,180,0,90,0,0,0,0\n";
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
      ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,1,0.005,0.001,0,0,0,0,0,0,0\n500GHz,1,0,180,1,1,0.005,-0.002,0.00001,0,0,0,0,0,0\n";
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
      Dataset sr = hfssfix::Mini500();
      sr.wg0 += "500GHz,0,0,180,3,4,0,0.001,0\n";   // a short row: 9 of 15 fields
      TempDir d3; ExpectG4Exception("BBR013", [&] { Load(d3, sr); }, "BBRHFSSData");
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
      TempDir d2; LoadChecked(d2, hfssfix::NineKey());                 // 3 far-field rows, 2 exit points per key
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
    {"sidecar_c5_outside_section", [] {
      const Dataset ds = hfssfix::WrapKeys("500", "0.01");   // exit points at Y = 10 mm
      auto p = hfssfix::SidecarFrom("c_500GHz", ds);
      p.yHalf = 5e-3;                                       // the declared section ends at 5 mm
      TempDir d; ExpectG4Exception("BBR025", [&] { LoadChecked(d, ds, &p); }, "BBRHFSSData",
                                   "5 exit point(s) lie outside the declared cross-section, e.g. (Y, Z) = (0.01, 1e-05) m");
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
  });
}
