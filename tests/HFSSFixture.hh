// HFSSFixture.hh — tiny synthetic HFSS datasets for the tests. Column orders
// match the producer's: waveguide.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,
//   Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag
// far_field.csv
//   Freq,Ephi,IWavePhi,IWaveTheta,Phi,Theta,rEphi_real,rEphi_imag,rEtheta_real,rEtheta_imag
#ifndef HFSSFixture_hh
#define HFSSFixture_hh
#include "BBRTestSupport.hh"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

// Three keys with distinct T0 (T1 = 0), one far-field row and one exit point each. Not a complete
// grid: its sidecar declares {0, 90} x {135, 180}, which lacks (90,135) (CDX-021).
inline Dataset ThreeKey() {
  Dataset d;
  d.ff0 = std::string(FF_HDR) + "500GHz,0,0,135,0,90,0,0,1,0\n500GHz,0,0,180,0,90,0,0,1,0\n500GHz,0,90,180,0,90,0,0,1,0\n";
  d.ff1 = std::string(FF_HDR) + "500GHz,1,0,135,0,90,1,0,0,0\n500GHz,1,0,180,0,90,1,0,0,0\n500GHz,1,90,180,0,90,1,0,0,0\n";
  d.wg0 = std::string(WG_HDR) + "500GHz,0,0,135,1,1,0,0,0,0,1,0,0,0,0\n500GHz,0,0,180,0.8,1,0,0,0,0,1,0,0,0,0\n500GHz,0,90,180,0.3,1,0,0,0,0,1,0,0,0,0\n";
  d.wg1 = std::string(WG_HDR) + "500GHz,1,0,135,0,1,0,0,0,0,0,1,0,0,0\n500GHz,1,0,180,0,1,0,0,0,0,0,1,0,0,0\n500GHz,1,90,180,0,1,0,0,0,0,0,1,0,0,0\n";
  return d;
}

// ThreeKey plus (90,135), T0 = 0.6: the complete {0, 90} x {135, 180} grid that ThreeKey's
// sidecar declares.
inline Dataset FourKey() {
  Dataset d = ThreeKey();
  d.ff0 += "500GHz,0,90,135,0,90,0,0,1,0\n";
  d.ff1 += "500GHz,1,90,135,0,90,1,0,0,0\n";
  d.wg0 += "500GHz,0,90,135,0.6,1,0,0,0,0,1,0,0,0,0\n";
  d.wg1 += "500GHz,1,90,135,0,1,0,0,0,0,0,1,0,0,0\n";
  return d;
}

// Six keys for the wrapper tests, at frequency <freq>GHz: {0, 45, 90} x {135, 180},
// the complete grid their sidecar declares (CDX-021). Polarization filters (T0 = 1,
// T1 = 0) at (0,135) and (0,180); T0 = T1 = 1 at the other four keys, where the
// orthogonal exit fields (y for Ephi=0, z for Ephi=1, rho = 0) make T = 1 for every
// polarization. Each key's far field is the full grid Phi {0, 30, 60} x Theta
// {30, 60, 90} (CDX-022), with intensity only at the direction that identifies the key:
//   (0,135) -> (Phi, Theta) = (30, 60)   (45,135) -> (0, 60)    (90,135) -> (60, 90)
//   (0,180) -> (0, 90)                   (45,180) -> (30, 90)   (90,180) -> (0, 30)
// One exit point per key at CSV (Y, Z) = (exitY m, 10 um).
inline Dataset WrapKeys(const std::string& freq = "500", const std::string& exitY = "0.001") {
  struct Row { const char* key; const char* phiTheta; const char* t0; const char* t1; };
  const Row rows[] = {
    {"0,135", "30,60", "1", "0"}, {"45,135", "0,60", "1", "1"},  {"90,135", "60,90", "1", "1"},
    {"0,180", "0,90", "1", "0"},  {"45,180", "30,90", "1", "1"}, {"90,180", "0,30", "1", "1"},
  };
  const std::string f = freq + "GHz,";
  Dataset d{FF_HDR, WG_HDR, FF_HDR, WG_HDR};
  for (const auto& r : rows) {
    const std::string k = std::string(r.key) + ",";
    for (const char* theta : {"30", "60", "90"})
      for (const char* phi : {"0", "30", "60"}) {
        const std::string pt = std::string(phi) + "," + theta;
        const bool lit = pt == r.phiTheta;
        d.ff0 += f + "0," + k + pt + (lit ? ",0,0,1,0\n" : ",0,0,0,0\n");
        d.ff1 += f + "1," + k + pt + (lit ? ",1,0,0,0\n" : ",0,0,0,0\n");
      }
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

// NineKey on a complete grid (CDX-022), for loads that must pass every check: per key the far
// field is Phi {-40, 10, 60} x Theta {40, 75, 110}, NineKey's three rows at (10, 40), (-40, 75)
// and (60, 110) with smaller amplitudes at the other six nodes, and the exit field the lattice
// Y {-2, 1} mm x Z {5, 10} um, NineKey's two points plus two without field. NineKey itself is
// not a complete grid: the C3/C4 tests that use it stop at the dataset-wide checks.
inline Dataset NineKeyGrid() {
  const char* keys[] = {"0,90", "45,90", "90,90", "0,135", "45,135", "90,135", "0,180", "45,180", "90,180"};
  Dataset d{FF_HDR, WG_HDR, FF_HDR, WG_HDR};
  for (const char* key : keys) {
    const std::string k = std::string(key) + ",";
    for (const char* theta : {"40", "75", "110"})
      for (const char* phi : {"-40", "10", "60"}) {
        const std::string pt = std::string(phi) + "," + theta;
        std::string a0 = "0.05,0.02,0.03,-0.01", a1 = "0.02,-0.03,0.04,0.01";
        if (pt == "10,40")  { a0 = "0.9,0.0,0.1,0.4";  a1 = "0.3,0.1,0.8,-0.2"; }
        if (pt == "-40,75") { a0 = "0.2,-0.3,0.7,0.1"; a1 = "-0.5,0.2,0.1,0.6"; }
        if (pt == "60,110") { a0 = "0.5,0.1,-0.6,0.3"; a1 = "0.4,0.4,-0.3,0.2"; }
        d.ff0 += "500GHz,0," + k + pt + "," + a0 + "\n";
        d.ff1 += "500GHz,1," + k + pt + "," + a1 + "\n";
      }
    d.wg0 += "500GHz,0," + k + "0.9,1,0,0.001,1e-05,0.7,0.1,0.4,0.2,0,0\n";
    d.wg0 += "500GHz,0," + k + "0.9,1,0,-0.002,5e-06,0.1,0.8,0.3,0,0.2,0.1\n";
    d.wg0 += "500GHz,0," + k + "0.9,1,0,0.001,5e-06,0,0,0,0,0,0\n";
    d.wg0 += "500GHz,0," + k + "0.9,1,0,-0.002,1e-05,0,0,0,0,0,0\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,0.001,1e-05,0.2,0.5,0.1,0,0.3,0\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,-0.002,5e-06,0.6,0.1,0.2,0.1,0,0.4\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,0.001,5e-06,0,0,0,0,0,0\n";
    d.wg1 += "500GHz,1," + k + "0.3,1,0,-0.002,1e-05,0,0,0,0,0,0\n";
  }
  return d;
}

namespace detail {
inline std::string Num17(double v) { std::ostringstream s; s << std::setprecision(17) << v; return s.str(); }

// One key (0,180), T0 = 1, T1 = 0, over the given exit lattice (Y, Z in m) with field only at the
// node (y_m, z_m) (y for Ephi=0, z for Ephi=1, as WrapKeys); a far field of one node (P_deg, T_deg)
// carrying all the intensity, so the direction and position draws are deterministic.
inline Dataset OnePoint(const std::vector<std::pair<double, double>>& lattice, double y_m, double z_m,
                        double P_deg, double T_deg) {
  if (std::find(lattice.begin(), lattice.end(), std::make_pair(y_m, z_m)) == lattice.end())
    throw std::invalid_argument("OnePoint: (" + Num17(y_m) + ", " + Num17(z_m) + ") is not a lattice node");
  const std::string key = "500GHz,", pt = Num17(P_deg) + "," + Num17(T_deg);
  Dataset d{FF_HDR, WG_HDR, FF_HDR, WG_HDR};
  d.ff0 += key + "0,0,180," + pt + ",0,0,1,0\n";
  d.ff1 += key + "1,0,180," + pt + ",1,0,0,0\n";
  for (const auto& [y, z] : lattice) {
    const bool lit = y == y_m && z == z_m;
    const std::string yz = Num17(y) + "," + Num17(z);
    d.wg0 += key + "0,0,180,1,1,0," + yz + (lit ? ",0,1,0,0,0,0\n" : ",0,0,0,0,0,0\n");
    d.wg1 += key + "1,0,180,0,1,0," + yz + (lit ? ",0,0,1,0,0,0\n" : ",0,0,0,0,0,0\n");
  }
  return d;
}
}  // namespace detail

// One exit point and one exit direction, for the wrapper's exit-start tests. The exit lattice is
// Y {-4.5, 0, 4.5} mm x Z {-25, 0, 25} um, the default declared section's walls and centre lines,
// so the section is unchanged and C2-C5 pass; (y_m, z_m) must be one of its nine nodes (a wall node
// is clamped 1 nm inside at load, OR-020). The far field is the single node (P_deg, T_deg). The
// frequency token is 500GHz; rename it with Relabel for another grid point.
inline Dataset OnePointOneDirection(double y_m, double z_m, double P_deg, double T_deg) {
  std::vector<std::pair<double, double>> lattice;
  for (double y : {-4.5e-3, 0., 4.5e-3})
    for (double z : {-2.5e-5, 0., 2.5e-5}) lattice.push_back({y, z});
  return detail::OnePoint(lattice, y_m, z_m, P_deg, T_deg);
}

// The same for a disc of radius R_m: the lattice {-R, 0, R}^2 inside the disc (5 nodes; the
// sidecar must declare the disc with outside_points "omitted").
inline Dataset OnePointOneDirectionDisc(double R_m, double y_m, double z_m, double P_deg, double T_deg) {
  const std::vector<std::pair<double, double>> lattice = {{-R_m, 0.}, {0., -R_m}, {0., 0.}, {0., R_m}, {R_m, 0.}};
  return detail::OnePoint(lattice, y_m, z_m, P_deg, T_deg);
}

// The dataset with its "500GHz," frequency tokens replaced by "<freq>GHz,".
inline Dataset Relabel(Dataset d, const std::string& freq) {
  for (std::string* s : {&d.ff0, &d.wg0, &d.ff1, &d.wg1})
    for (std::size_t p = 0; (p = s->find("500GHz,", p)) != std::string::npos; p += freq.size() + 4)
      s->replace(p, 7, freq + "GHz,");
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
  double lengthMm = 4.0;   // geometry.extent_mm.p; l and g follow from the section
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
    if ((v.size() != 15 && v.size() != 16) || !detail::Num(v[2], ip) || !detail::Num(v[3], it) ||
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
  // extent_mm: p the crack length, l and g the section's full widths (a disc: its diameter), in mm,
  // as the producer writes them; the parse checks a rectangle's l and g against its section.
  const double extL = p.disc ? 2e3 * p.radius : 2e3 * p.yHalf, extG = p.disc ? 2e3 * p.radius : 2e3 * p.zHalf;
  const std::string extent = "{\"p\": " + Fmt(p.lengthMm) + ", \"l\": " + Fmt(extL) + ", \"g\": " + Fmt(extG) + "}";
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
    << "    \"exit_cs\": {\"name\": \"outgoing_cs\", \"origin\": \"exit_face_center\", \"origin_mm_global\": [0, 0, 4], \"x\": [0, 0, 1], \"y\": [0, 1, 0], \"z\": [-1, 0, 0]},\n"
    << "    \"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, 1, 0], \"z\": [0, 0, 1]}},\n"
    << "  \"excitation\": {\"coordinate_system\": \"global\", \"incidence_convention\": \"arrival_direction\",\n"
    << "    \"normal_entry_theta_deg\": 180, \"incident_phi_deg\": " << list(p.phis)
    << ", \"incident_theta_deg\": " << list(p.thetas) << ",\n"
    << "    \"polarization_convention\": \"Ephi=0: E_theta=1; Ephi=1: E_phi=1\", \"origin_mm_global\": [0, 0, 0]},\n"
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
    << "  \"geometry\": {\"shape\": \"" << (p.disc ? "cylinder" : "box") << "\", \"extent_mm\": " << extent << ",\n"
    << "    \"bounding_box_mm\": [-0.026, -5, 0, 0.026, 5, 4]},\n"
    << "  \"modes\": {\"cutoff_ghz\": " << Fmt(cut) << ", \"mode\": \"" << mode << "\", \"polarization_filter_limit_ghz\": "
    << Fmt(filter) << ", \"propagating_count\": 0}\n"
    << "}\n";
  return j.str();
}

// --- schema 2.0 -----------------------------------------------------------------
// Sixteen-column waveguide rows (...,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag,InterceptFactor)
// and a 2.0 sidecar carrying every field of the 2.0 field list. SidecarJson2 starts from the 1.0
// text of the same SidecarParams and rewrites what 2.0 changes, so the shared blocks cannot drift.
const char* const WG_HDR2 = "Freq,Ephi,IWavePhi,IWaveTheta,OutgoingPower,IngoingPower,X,Y,Z,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag,InterceptFactor\n";

// tflux holds one {t0, t1, dark0, dark1} per declared key in declared order (phi-major, then theta;
// dark as 0/1); empty means t0 = 0.6, t1 = 0.3, nothing dark. yStep/zStep of 0 are derived from the
// declared lattice. emitter and floor feed transmittance.cell_90_assumes.
struct Schema2Params {
  SidecarParams base;
  std::vector<std::array<double, 4>> tflux;
  double yStep = 0., zStep = 0.;
  std::string emitter = "uniform_theta";
  double floor = 0.01;
};

namespace detail {
constexpr double kPiF = 3.14159265358979323846;
// A fixture-local copy of the HF-049 rule, not shared with the library: the reference for the
// parser's own copy. Delta = (smallest declared theta above 90) - 90; -1 if there is none or no 90.
inline double DeltaF(const std::vector<double>& th) {
  bool has90 = false;
  double above = 1e300;
  for (double t : th) {
    if (std::abs(t - 90.) <= 1e-9) has90 = true;
    else if (t > 90.) above = std::min(above, t);
  }
  return (has90 && above < 1e300) ? above - 90. : -1.;
}
inline double Cell90F(const std::vector<double>& th) {
  const double d = DeltaF(th);
  if (d < 0.) return -1.;
  const double a0 = kPiF / 2 - d / 2 * kPiF / 180., af = std::acos(0.01), L = kPiF / 2 - a0;
  const auto prim = [](double a) { return std::log(1. / std::cos(a) + std::tan(a)); };
  const double I = af > a0 ? (prim(af) - prim(a0)) + (kPiF / 2 - af) / 0.01 : (kPiF / 2 - a0) / 0.01;
  return L / I;
}
inline double InterceptF(double theta, const std::vector<double>& th) {
  return std::abs(theta - 90.) <= 1e-9 ? Cell90F(th) : std::abs(std::cos(theta * kPiF / 180.));
}
}  // namespace detail

inline std::string SidecarJson2(const Schema2Params& p) {
  using nlohmann::json;
  constexpr double c = 299792458.0;
  const SidecarParams& b = p.base;
  json j = json::parse(SidecarJson(b));
  j["schema_version"] = "2.0";
  j["model"] = {{"kind", "screen_free_abc"}};
  j["conventions"] = {{"time_dependence", "exp(+j omega t)"}, {"phasor", "peak"}};

  const double cell90 = detail::Cell90F(b.thetas), delta = detail::DeltaF(b.thetas);
  json cKey = json::array(), perKey = json::array(), scalePer = json::array(), ffRatio = json::array();
  std::size_t k = 0;
  for (double phi : b.phis)
    for (double theta : b.thetas) {
      const std::array<double, 4> t = k < p.tflux.size() ? p.tflux[k] : std::array<double, 4>{0.6, 0.3, 0., 0.};
      const bool d0 = t[2] != 0., d1 = t[3] != 0.;
      const double ck = detail::InterceptF(theta, b.thetas);
      cKey.push_back({{"phi_deg", phi}, {"theta_deg", theta}, {"c_k", ck}});
      perKey.push_back({{"phi_deg", phi}, {"theta_deg", theta}, {"c_k", ck},
                        {"polarizations", {{"0", {{"t_flux", d0 ? 0. : t[0]}, {"outgoing_power_abs_re_s_w", 1.05e-9 * t[0]}}},
                                           {"1", {{"t_flux", d1 ? 0. : t[1]}, {"outgoing_power_abs_re_s_w", 1.05e-9 * t[1]}}}}}});
      scalePer.push_back({{"phi_deg", phi}, {"theta_deg", theta}, {"w0", d0 ? 1. : 0.9}, {"w1", d1 ? 1. : 0.8},
                          {"dark0", d0}, {"dark1", d1}});
      ffRatio.push_back({{"phi_deg", phi}, {"theta_deg", theta},
                         {"polarizations", {{"0", d0 ? json(nullptr) : json(0.179)}, {"1", d1 ? json(nullptr) : json(0.179)}}}});
      ++k;
    }
  json band = json::array({90. - delta / 2, 90.});
  j["transmittance"] = {
      {"definition", "OutgoingPower / (IngoingPower * InterceptFactor)"},
      {"numerator", "exit_flux_re_s_dot_n"},
      {"incoming", "area_times_intensity_undivided"},
      {"above_one", "per_photon_cap"},
      {"c_key", cKey},
      {"cell_90", cell90},
      {"cell_90_assumes", {{"emitter", p.emitter}, {"cos_floor", p.floor}, {"band_deg", band}}},
      {"per_key", perKey}};

  json& xf = j["exit_field"];
  const double ys = p.yStep > 0. ? p.yStep : (b.yCount > 1 ? (b.yMax - b.yMin) / (b.yCount - 1) : 0.);
  const double zs = p.zStep > 0. ? p.zStep : (b.zCount > 1 ? (b.zMax - b.zMin) / (b.zCount - 1) : 0.);
  xf["grid"]["y_e"]["step"] = ys;
  xf["grid"]["z_e"]["step"] = zs;
  xf["columns"] = {"Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower", "X", "Y", "Z", "Ex_real",
                   "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor"};
  xf["tangential_components"] = {"Ex", "Ey"};
  xf["quadrature"] = b.disc ? "uniform_disc" : "trapezoid";
  xf["sampler"] = "abs_e_t_squared_times_quadrature";
  xf["scale"] = {{"rule", "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w"}, {"per_key", scalePer}};
  xf["column_order_check"] = {{"rule", "reference_exporter_interleaved_v1"}, {"method", "reference"},
                              {"reference", {{"dataset_id", "InfParallelPlate_crack1Rohan"}, {"frequency_ghz", 500},
                                             {"passed", true}}}};
  xf["plane_tolerance_m"] = 1e-12;

  j["far_field"]["normalization"] = "pattern_only";
  j["far_field"]["ff_over_flux"] = ffRatio;

  json& md = j["modes"];
  md.erase("propagating_count");
  const double filter = md.at("polarization_filter_limit_ghz").get<double>();
  json regime = json::array();
  if (b.disc) {
    for (const auto& [name, zero] : std::vector<std::pair<const char*, double>>{
             {"TE11", 1.8411837813}, {"TM01", 2.404825557695773}, {"TE21", 3.054236928227140}})
      regime.push_back({{"mode", name}, {"cutoff_ghz", zero * c / (2 * detail::kPiF * b.radius) / 1e9}});
  } else {
    const double a = 2 * b.yHalf, bb = 2 * b.zHalf;
    regime.push_back({{"mode", "TE10"}, {"cutoff_ghz", c / (2 * a) / 1e9}});
    for (int n = 1; n * c / (2 * bb) / 1e9 <= 20000. + 1e-9; ++n)
      regime.push_back({{"mode", n < 10 ? "TE0" + std::to_string(n) : "TE0," + std::to_string(n)},
                        {"cutoff_ghz", n * c / (2 * bb) / 1e9}});
  }
  // Ascending by cutoff, as the producer's ordered list (differs from "TE10 first" only if a < b).
  std::stable_sort(regime.begin(), regime.end(), [](const json& x, const json& y) {
    return x.at("cutoff_ghz").get<double>() < y.at("cutoff_ghz").get<double>();
  });
  md["regime_cutoffs_ghz"] = regime;

  j["cross_term"] = {{"density", "e_t_quadrature"},
                     {"definition", "Re rho = sum_i w_i Re(E_t0,i . conj E_t1,i) / sqrt(sum_i w_i abs(E_t0,i)^2 * sum_i w_i abs(E_t1,i)^2)"},
                     {"valid_below_ghz", filter},
                     {"above_limit", "warn_once"}};
  j["transmittance"].erase("incoming_includes_cos_theta");
  return j.dump(2) + "\n";
}

// The nine keys of NineKeyGrid with 16-column rows: the same far field, a 3 x 3 exit lattice
// y {-1e-3, 0, 1e-3} m x z {-2.5e-5, 0, 2.5e-5} m, IngoingPower 1e-9, OutgoingPower = t_p * 1e-9
// with t0 = 0.6 and t1 = 0.3 except Ephi=1 at (0, 180), which is dark (OutgoingPower 0 with
// nonzero stored fields), InterceptFactor = the key's c_k, and the six field values of point i
// equal to (i+1 ... i+6) * 1e-3 for Ephi=0 and in reverse order for Ephi=1.
inline Dataset NineKeyGrid2() {
  const Dataset g = NineKeyGrid();
  const std::vector<double> th{90., 135., 180.};
  Dataset d{g.ff0, WG_HDR2, g.ff1, WG_HDR2};
  const double ys[] = {-1e-3, 0., 1e-3}, zs[] = {-2.5e-5, 0., 2.5e-5};
  for (int ph : {0, 45, 90})
    for (int t : {90, 135, 180}) {
      const std::string key = detail::Fmt(ph) + "," + detail::Fmt(t) + ",";
      const std::string ck = detail::Fmt(detail::InterceptF(t, th));
      for (int e : {0, 1}) {
        const bool dark = e == 1 && ph == 0 && t == 180;
        const double tp = e == 0 ? 0.6 : 0.3;
        std::string& wg = e == 0 ? d.wg0 : d.wg1;
        int i = 0;
        for (double y : ys)
          for (double z : zs) {
            wg += "500GHz," + std::to_string(e) + "," + key + detail::Fmt(dark ? 0. : tp * 1e-9) + ",1e-09,0," +
                  detail::Fmt(y) + "," + detail::Fmt(z);
            for (int f = 0; f < 6; ++f) wg += "," + detail::Fmt((e == 0 ? i + 1 + f : i + 6 - f) * 1e-3);
            wg += "," + ck + "\n";
            ++i;
          }
      }
    }
  return d;
}

inline Schema2Params Schema2From(const std::string& stem, const Dataset& d) {
  Schema2Params p;
  p.base = SidecarFrom(stem, d);
  return p;
}

// The params of NineKeyGrid2: (0, 180) Ephi=1 is dark.
inline Schema2Params NineKeyParams2(const std::string& stem = "crack_500GHz") {
  Schema2Params p = Schema2From(stem, NineKeyGrid2());
  p.tflux.assign(p.base.phis.size() * p.base.thetas.size(), std::array<double, 4>{0.6, 0.3, 0., 0.});
  p.tflux[2] = {0.6, 0., 0., 1.};   // phi 0: theta 90, 135, 180
  return p;
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

// Writes the CSVs of d and the 2.0 sidecar of p.
inline std::filesystem::path WriteDataset2(const std::filesystem::path& root, const std::string& stem,
                                           const Dataset& d, const Schema2Params& p) {
  const auto wg = WriteDataset(root, stem, d, false);
  bbrtest::WriteFile(wg / (stem + ".dataset.json"), SidecarJson2(p));
  return wg;
}

}  // namespace hfssfix
#endif
