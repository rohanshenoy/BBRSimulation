// HFSSFixture.hh — tiny synthetic HFSS datasets for the tests. Column orders
// match the producer's: waveguide.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,
//   Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag
// far_field.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag
#ifndef HFSSFixture_hh
#define HFSSFixture_hh
#include "BBRTestSupport.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
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

// --- sidecars (schema 1.0) ---------------------------------------------------
// BBRCrackLibrary refuses a dataset without <stem>.dataset.json (BBR024), so
// WriteDataset writes one, derived from the dataset's own rows by SidecarFrom.
// The default exit section, half-extents 4.5 mm x 25 um, sits inside the
// 5 mm x 26 um crack of WrapperWorld and testPhysicsList (the F12 fit check)
// and grows to cover rows that reach further (WrapKeys("500", "0.01")). The
// modes follow from the section as BBRDatasetSidecar::Parse requires (F13).
struct SidecarParams {
  std::string id, label;
  double freqGHz = 0.;
  std::vector<double> phis{0.}, thetas{180.};
  double thMin = 90., thMax = 90.; int thCount = 1;
  double phMin = 0., phMax = 0.; int phCount = 1; int ffPerKey = 1;
  double yMin = 0., yMax = 0.; int yCount = 1;
  double zMin = 0., zMax = 0.; int zCount = 1; int exitPerKey = 1;
  bool disc = false;
  double yHalf = 4.5e-3, zHalf = 2.5e-5, radius = 0.;
  std::string outsidePoints = "none";
};

namespace detail {
inline std::vector<std::vector<std::string>> Rows(const std::string& csv) {
  std::vector<std::vector<std::string>> out;
  std::istringstream in(csv);
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (header) { header = false; continue; }
    if (line.empty()) continue;
    std::vector<std::string> v;
    std::stringstream ss(line);
    std::string t;
    while (std::getline(ss, t, ',')) v.push_back(t);
    out.push_back(v);
  }
  return out;
}
inline bool Num(const std::string& s, double& v) {
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return !s.empty() && end && *end == '\0' && std::isfinite(v);
}
inline std::string Fmt(double v) { std::ostringstream s; s << std::setprecision(17) << v; return s.str(); }
// "<id>_<number>GHz": the stems BBRCrackLibrary discovers.
inline bool FrequencyStem(const std::string& stem) {
  const auto us = stem.rfind('_');
  if (us == std::string::npos || us == 0 || stem.size() < us + 5) return false;
  const std::string label = stem.substr(us + 1);
  char* end = nullptr;
  const double f = std::strtod(label.c_str(), &end);
  return f > 0. && end && std::string(end) == "GHz";
}
}  // namespace detail

// Rows that do not parse are skipped, so deliberately broken fixtures still get one.
inline SidecarParams SidecarFrom(const std::string& stem, const Dataset& d) {
  SidecarParams p;
  const auto us = stem.rfind('_');
  p.id = stem.substr(0, us);
  p.label = stem.substr(us + 1);
  p.freqGHz = std::strtod(p.label.c_str(), nullptr);
  std::set<double> phis, thetas, P, T, Y, Z;
  std::map<std::pair<double, double>, int> ffN, wgN;
  for (const auto& v : detail::Rows(d.ff0.empty() ? d.ff1 : d.ff0)) {
    double ip, it, ph, th;
    if (v.size() != 10 || !detail::Num(v[2], ip) || !detail::Num(v[3], it) ||
        !detail::Num(v[4], ph) || !detail::Num(v[5], th)) continue;
    phis.insert(ip); thetas.insert(it); P.insert(ph); T.insert(th); ++ffN[{ip, it}];
  }
  for (const auto& v : detail::Rows(d.wg0.empty() ? d.wg1 : d.wg0)) {
    double ip, it, y, z;
    if (v.size() != 15 || !detail::Num(v[2], ip) || !detail::Num(v[3], it) ||
        !detail::Num(v[7], y) || !detail::Num(v[8], z)) continue;
    phis.insert(ip); thetas.insert(it); Y.insert(y); Z.insert(z); ++wgN[{ip, it}];
  }
  if (!phis.empty()) { p.phis.assign(phis.begin(), phis.end()); p.thetas.assign(thetas.begin(), thetas.end()); }
  if (!P.empty()) {
    p.phMin = *P.begin(); p.phMax = *P.rbegin(); p.phCount = int(P.size());
    p.thMin = *T.begin(); p.thMax = *T.rbegin(); p.thCount = int(T.size());
    p.ffPerKey = ffN.begin()->second;
  }
  if (!Y.empty()) {
    p.yMin = *Y.begin(); p.yMax = *Y.rbegin(); p.yCount = int(Y.size());
    p.zMin = *Z.begin(); p.zMax = *Z.rbegin(); p.zCount = int(Z.size());
    p.exitPerKey = wgN.begin()->second;
    p.yHalf = std::max(p.yHalf, std::max(std::abs(p.yMin), std::abs(p.yMax)));
    p.zHalf = std::max(p.zHalf, std::max(std::abs(p.zMin), std::abs(p.zMax)));
  }
  return p;
}

inline std::string SidecarJson(const SidecarParams& p) {
  constexpr double c = 299792458.0, pi = 3.14159265358979323846;
  std::string mode;
  double cut, filter;
  if (p.disc) {
    mode = "TE11"; cut = 1.8411837813 * c / (2 * pi * p.radius) / 1e9; filter = cut;
  } else {
    const double a = 2 * p.yHalf, b = 2 * p.zHalf;
    if (a >= b) { mode = "TE10"; cut = c / (2 * a) / 1e9; } else { mode = "TE01"; cut = c / (2 * b) / 1e9; }
    filter = c / (2 * b) / 1e9;
  }
  auto list = [](const std::vector<double>& v) {
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + detail::Fmt(v[i]);
    return s + "]";
  };
  using detail::Fmt;
  const std::string section = p.disc
    ? "{\"shape\": \"disc\", \"radius_m\": " + Fmt(p.radius) + "}"
    : "{\"shape\": \"rectangle\", \"y_e_half_m\": " + Fmt(p.yHalf) + ", \"z_e_half_m\": " + Fmt(p.zHalf) + "}";
  std::ostringstream j;
  j << "{\n"
    << "  \"schema_version\": \"1.0\",\n"
    << "  \"dataset_id\": \"" << p.id << "\",\n"
    << "  \"frequency_ghz\": " << Fmt(p.freqGHz) << ",\n"
    << "  \"frequency_label\": \"" << p.label << "\",\n"
    << "  \"provenance\": {\"producer\": \"tests/HFSSFixture.hh\"},\n"
    << "  \"frames\": {\"canonical\": \"p,l,g\", \"pose_rule\": \"canonical-z\",\n"
    << "    \"hfss_global_axes_in_canonical\": {\"x\": [0, 0, -1], \"y\": [0, 1, 0], \"z\": [1, 0, 0]},\n"
    << "    \"entrance_face\": {\"axis\": \"z\", \"side\": \"min\"}, \"exit_face\": {\"axis\": \"z\", \"side\": \"max\"},\n"
    << "    \"entrance_outward_normal_global\": [0, 0, -1], \"exit_outward_normal_global\": [0, 0, 1],\n"
    << "    \"exit_cs\": {\"name\": \"outgoing_cs\", \"origin\": \"exit_face_center\", \"x\": [0, 0, 1], \"y\": [0, 1, 0], \"z\": [-1, 0, 0]},\n"
    << "    \"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, 1, 0], \"z\": [0, 0, 1]}},\n"
    << "  \"excitation\": {\"coordinate_system\": \"global\", \"incidence_convention\": \"arrival_direction\",\n"
    << "    \"normal_entry_theta_deg\": 180, \"incident_phi_deg\": " << list(p.phis)
    << ", \"incident_theta_deg\": " << list(p.thetas) << ",\n"
    << "    \"polarization_convention\": \"Ephi=0: E_theta=1; Ephi=1: E_phi=1\"},\n"
    << "  \"far_field\": {\"coordinate_system\": \"exit_cs\", \"definition\": \"Theta-Phi\", \"component_basis\": \"spherical_in_exit_cs\",\n"
    << "    \"theta_deg\": {\"min\": " << Fmt(p.thMin) << ", \"max\": " << Fmt(p.thMax) << ", \"count\": " << p.thCount << "},\n"
    << "    \"phi_deg\": {\"min\": " << Fmt(p.phMin) << ", \"max\": " << Fmt(p.phMax) << ", \"count\": " << p.phCount << "},\n"
    << "    \"points_per_key\": " << p.ffPerKey << ",\n"
    << "    \"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"Phi\", \"Theta\", \"rEphi_real\", \"rEphi_imag\", \"rEtheta_real\", \"rEtheta_imag\"]},\n"
    << "  \"exit_field\": {\"coordinate_system\": \"exit_cs\", \"points_in_si\": true, \"field_in_ref_cs\": false,\n"
    << "    \"field_components_frame\": \"hfss_global\", \"plane\": \"x_e=0\",\n"
    << "    \"grid\": {\"y_e\": {\"min\": " << Fmt(p.yMin) << ", \"max\": " << Fmt(p.yMax) << ", \"count\": " << p.yCount << "}, "
    << "\"z_e\": {\"min\": " << Fmt(p.zMin) << ", \"max\": " << Fmt(p.zMax) << ", \"count\": " << p.zCount << "}},\n"
    << "    \"cross_section\": " << section << ",\n"
    << "    \"outside_points\": \"" << p.outsidePoints << "\", \"rim_points\": \"included\", \"points_per_key_retained\": " << p.exitPerKey << ",\n"
    << "    \"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"OutgoingPower\", \"IngoingPower\", \"X\", \"Y\", \"Z\", \"Ex_real\", \"Ey_real\", \"Ez_real\", \"Ex_imag\", \"Ey_imag\", \"Ez_imag\"]},\n"
    << "  \"transmittance\": {\"definition\": \"outgoing_power_w / incoming_power_w, first row per key\", \"incoming_includes_cos_theta\": false},\n"
    << "  \"symmetry\": {\"mirror_l\": true, \"mirror_g\": true, \"end_to_end\": true, \"rotational\": " << (p.disc ? "true" : "false") << "},\n"
    << "  \"boundaries\": {\"entrance\": \"radiation\", \"exit\": \"radiation\", \"walls\": \"PEC\"},\n"
    << "  \"geometry\": {\"shape\": \"" << (p.disc ? "cylinder" : "box") << "\", \"extent_mm\": {\"p\": 4, \"l\": 10, \"g\": 0.052}},\n"
    << "  \"modes\": {\"cutoff_ghz\": " << Fmt(cut) << ", \"mode\": \"" << mode << "\", \"polarization_filter_limit_ghz\": "
    << Fmt(filter) << ", \"propagating_count\": 0}\n"
    << "}\n";
  return j.str();
}

// Writes <root>/waveguides/<stem>_Ephi={0,1}/{far_field,waveguide}.csv and, for
// a <id>_<freq>GHz stem with withSidecar, <stem>.dataset.json from SidecarFrom.
// An empty string skips that file (to provoke BBR001/BBR002).
inline std::filesystem::path WriteDataset(const std::filesystem::path& root,
                                          const std::string& stem, const Dataset& d,
                                          bool withSidecar = true) {
  const auto wg = root / "waveguides";
  if (!d.ff0.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=0") / "far_field.csv", d.ff0);
  if (!d.wg0.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=0") / "waveguide.csv", d.wg0);
  if (!d.ff1.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=1") / "far_field.csv", d.ff1);
  if (!d.wg1.empty()) bbrtest::WriteFile(wg / (stem + "_Ephi=1") / "waveguide.csv", d.wg1);
  std::filesystem::create_directories(wg / (stem + "_Ephi=0"));
  std::filesystem::create_directories(wg / (stem + "_Ephi=1"));
  if (withSidecar && detail::FrequencyStem(stem))
    bbrtest::WriteFile(wg / (stem + ".dataset.json"), SidecarJson(SidecarFrom(stem, d)));
  return wg;
}

}  // namespace hfssfix
#endif
