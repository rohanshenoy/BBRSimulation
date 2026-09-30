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
