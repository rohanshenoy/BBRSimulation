// HFSSFixture.hh — tiny synthetic HFSS datasets for the tests. Column orders
// match the producer's: waveguide.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,
//   Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag
// far_field.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag
#ifndef HFSSFixture_hh
#define HFSSFixture_hh
#include "BBRTestSupport.hh"
#include <filesystem>
#include <string>

namespace hfssfix {

struct Dataset { std::string ff0, wg0, ff1, wg1; };

const char* const WG_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag\n";
const char* const FF_HDR = "Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag\n";

// One key (0,180); T0 = 3/4, T1 = 1/4; one far-field row straight out
// (Theta 90, Phi 0); one exit point at the face centre, field along y (Ephi=0)
// and z (Ephi=1) so the two exit fields are orthogonal (rho = 0).
inline Dataset Mini500() {
  Dataset d;
  d.ff0 = std::string(FF_HDR) + "500GHz,0,0,180,0,90,0,0,1,0\n";
  d.ff1 = std::string(FF_HDR) + "500GHz,1,0,180,0,90,1,0,0,0\n";
  d.wg0 = std::string(WG_HDR) + "500GHz,0,0,180,3,4,0,0,0,0,1,0,0,0,0\n";
  d.wg1 = std::string(WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,0,0,1,0,0,0\n";
  return d;
}

// Three keys with distinct T0 (T1 = 0), one far-field row and one exit point each.
inline Dataset ThreeKey() {
  Dataset d;
  d.ff0 = std::string(FF_HDR) + "500GHz,0,0,135,0,90,0,0,1,0\n500GHz,0,0,180,0,90,0,0,1,0\n500GHz,0,90,180,0,90,0,0,1,0\n";
  d.ff1 = std::string(FF_HDR) + "500GHz,1,0,135,0,90,1,0,0,0\n500GHz,1,0,180,0,90,1,0,0,0\n500GHz,1,90,180,0,90,1,0,0,0\n";
  d.wg0 = std::string(WG_HDR) + "500GHz,0,0,135,1,1,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,0.8,1,0,0,0,0,1,0,0,0,0\n500GHz,0,90,180,0.3,1,0,0,0,0,1,0,0,0,0\n";
  d.wg1 = std::string(WG_HDR) + "500GHz,1,0,135,0,1,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,0,1,0,0,0,0,0,1,0,0,0\n500GHz,1,90,180,0,1,0,0,0,0,0,1,0,0,0\n";
  return d;
}

// Five keys for the wrapper tests, at frequency <freq>GHz: the three keys of
// ThreeKey plus (90,135) and (45,135). Polarization filters (T0 = 1, T1 = 0) at
// (0,135) and (0,180); T0 = T1 = 1 at (45,135), (90,135) and (90,180), where the
// orthogonal exit fields (y for Ephi=0, z for Ephi=1, rho = 0) make T = 1 for
// every polarization. One far-field row per key, at a direction that identifies
// the key:
//   (0,135) -> (Phi, Theta) = (30, 60)   (45,135) -> (0, 60)   (90,135) -> (60, 90)
//   (0,180) -> (0, 90)                   (90,180) -> (0, 45)
// One exit point per key at CSV (Y, Z) = (exitY m, 10 um).
inline Dataset WrapKeys(const std::string& freq = "500", const std::string& exitY = "0.001") {
  struct Row { const char* key; const char* phiTheta; const char* t0; const char* t1; };
  const Row rows[] = {
    {"0,135", "30,60", "1", "0"}, {"45,135", "0,60", "1", "1"}, {"90,135", "60,90", "1", "1"},
    {"0,180", "0,90", "1", "0"},  {"90,180", "0,45", "1", "1"},
  };
  const std::string f = freq + "GHz,";
  Dataset d{FF_HDR, WG_HDR, FF_HDR, WG_HDR};
  for (const auto& r : rows) {
    const std::string k = std::string(r.key) + ",", pt = std::string(r.phiTheta) + ",";
    d.ff0 += f + "0," + k + pt + "0,0,1,0\n";
    d.ff1 += f + "1," + k + pt + "1,0,0,0\n";
    d.wg0 += f + "0," + k + r.t0 + ",1,0," + exitY + ",1e-05,0,1,0,0,0,0\n";
    d.wg1 += f + "1," + k + r.t1 + ",1,0," + exitY + ",1e-05,0,0,1,0,0,0\n";
  }
  return d;
}

// Nine keys (IWavePhi 0/45/90 x IWaveTheta 90/135/180), every key T0 = 0.9,
// T1 = 0.3, three far-field rows with complex amplitudes that differ between
// Ephi=0 and Ephi=1, and two exit points with different complex fields (so
// rho != 0 and the P1 cross term is active). For the mirror-equivariance test.
inline Dataset NineKey() {
  const char* keys[] = {"0,90", "45,90", "90,90", "0,135", "45,135", "90,135", "0,180", "45,180", "90,180"};
  Dataset d{FF_HDR, WG_HDR, FF_HDR, WG_HDR};
  for (const char* key : keys) {
    const std::string k = std::string(key) + ",";
    d.ff0 += "500GHz,0," + k + "20,40,0.9,0.0,0.1,0.4\n";
    d.ff0 += "500GHz,0," + k + "-40,70,0.2,-0.3,0.7,0.1\n";
    d.ff0 += "500GHz,0," + k + "60,110,0.5,0.1,-0.6,0.3\n";
    d.ff1 += "500GHz,1," + k + "20,40,0.3,0.1,0.8,-0.2\n";
    d.ff1 += "500GHz,1," + k + "-40,70,-0.5,0.2,0.1,0.6\n";
    d.ff1 += "500GHz,1," + k + "60,110,0.4,0.4,-0.3,0.2\n";
    d.wg0 += "500GHz,0," + k + "0.9,1,0,0.001,1e-05,0.7,0.1,0.4,0.2,0,0\n";
    d.wg0 += "500GHz,0," + k + "0.9,1,0,-0.002,5e-06,0.1,0.8,0.3,0,0.2,0.1\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,0.001,1e-05,0.2,0.5,0.1,0,0.3,0\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,-0.002,5e-06,0.6,0.1,0.2,0.1,0,0.4\n";
  }
  return d;
}

// Writes <root>/waveguides/<stem>_Ephi={0,1}/{far_field,waveguide}.csv.
// An empty string skips that file (to provoke BBR001/BBR002).
inline std::filesystem::path WriteDataset(const std::filesystem::path& root,
                                          const std::string& stem, const Dataset& d) {
  const auto wg = root / "waveguides";
  if (!d.ff0.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=0") / "far_field.csv", d.ff0);
  if (!d.wg0.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=0") / "waveguide.csv", d.wg0);
  if (!d.ff1.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=1") / "far_field.csv", d.ff1);
  if (!d.wg1.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=1") / "waveguide.csv", d.wg1);
  std::filesystem::create_directories(wg / (stem + "_Ephi=0"));
  std::filesystem::create_directories(wg / (stem + "_Ephi=1"));
  return wg;
}

}  // namespace hfssfix
#endif
