#include "BBRDatasetSidecar.hh"

#include "BBRHFSSData.hh"

#include "G4Exception.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "geomdefs.hh"

#include "nlohmann/json.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

using json = nlohmann::json;

namespace {

using Vec = std::array<double, 3>;

constexpr double kC = 299792458.0;           // m/s, exact
constexpr double kX11Prime = 1.8411837813;   // first zero of J1': TE11 of a round guide
constexpr double kPi = 3.14159265358979323846;
constexpr double kAxisTol = 1e-9;
constexpr double kModeRelTol = 1e-6;

const char* const kPolarization = "Ephi=0: E_theta=1; Ephi=1: E_phi=1";
const char* const kTransmittance = "outgoing_power_w / incoming_power_w, first row per key";
const std::vector<std::string> kFarFieldColumns = {"Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta",
                                                   "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"};
const std::vector<std::string> kExitColumns = {"Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower",
                                               "IngoingPower", "X", "Y", "Z", "Ex_real", "Ey_real", "Ez_real",
                                               "Ex_imag", "Ey_imag", "Ez_imag"};
// The frame the sampler implements (validation/README.md, Dataset sidecars).
const Vec kHfssX{0, 0, -1}, kHfssY{0, 1, 0}, kHfssZ{1, 0, 0};
const Vec kP{1, 0, 0}, kL{0, 1, 0}, kG{0, 0, 1};

[[noreturn]] void Fail(const char* code, const std::string& path, const std::string& what) {
  G4ExceptionDescription ed;
  ed << "HFSS dataset sidecar " << path << ": " << what;
  G4Exception("BBRDatasetSidecar", code, FatalException, ed);
  throw std::logic_error("a fatal G4Exception returned");   // only under a non-aborting handler
}

std::string Fmt(double v) { std::ostringstream s; s.precision(12); s << v; return s.str(); }
std::string Fmt(const Vec& v) { return "[" + Fmt(v[0]) + ", " + Fmt(v[1]) + ", " + Fmt(v[2]) + "]"; }
double Dot(const Vec& a, const Vec& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec Cross(const Vec& a, const Vec& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
bool Near(const Vec& a, const Vec& b) {
  for (int i = 0; i < 3; ++i) if (std::abs(a[i] - b[i]) > kAxisTol) return false;
  return true;
}
bool RelNear(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// Field access with type checks: a missing or mistyped field is BBR024.
struct Reader {
  const std::string& path;
  const json& At(const json& j, const std::string& key, const std::string& where) const {
    if (!j.is_object() || !j.contains(key)) Fail("BBR024", path, "missing field " + where + key);
    return j.at(key);
  }
  std::string Str(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    if (!v.is_string()) Fail("BBR024", path, where + key + " must be a string");
    return v.get<std::string>();
  }
  double Num(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    if (!v.is_number() || !std::isfinite(v.get<double>())) Fail("BBR024", path, where + key + " must be a finite number");
    return v.get<double>();
  }
  int Int(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    // Range-checked before narrowing: a count beyond int is malformed, not wrapped.
    const bool inRange = v.is_number_unsigned() ? v.get<std::uint64_t>() <= std::uint64_t(INT_MAX)
                       : v.is_number_integer() ? (v.get<std::int64_t>() >= INT_MIN && v.get<std::int64_t>() <= INT_MAX)
                                               : false;
    if (!inRange) Fail("BBR024", path, where + key + " must be an integer within the int range");
    return static_cast<int>(v.get<std::int64_t>());
  }
  bool Bool(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    if (!v.is_boolean()) Fail("BBR024", path, where + key + " must be true or false");
    return v.get<bool>();
  }
  Vec Vec3(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    if (!v.is_array() || v.size() != 3) Fail("BBR024", path, where + key + " must be a 3-vector");
    Vec out{};
    for (std::size_t i = 0; i < 3; ++i) {
      if (!v[i].is_number()) Fail("BBR024", path, where + key + " must be a 3-vector of numbers");
      out[i] = v[i].get<double>();
    }
    return out;
  }
  std::vector<double> Nums(const json& j, const std::string& key, const std::string& where) const {
    const json& v = At(j, key, where);
    if (!v.is_array() || v.empty()) Fail("BBR024", path, where + key + " must be a non-empty array of numbers");
    std::vector<double> out;
    for (const auto& x : v) {
      if (!x.is_number()) Fail("BBR024", path, where + key + " must be a non-empty array of numbers");
      out.push_back(x.get<double>());
    }
    return out;
  }
  // A CSV column list: not an array is a mistyped field (BBR024); an array that is
  // not exactly `want`, non-string entries included, is BBR025, as in the Python check.
  std::vector<std::string> Columns(const json& j, const std::string& key, const std::string& where,
                                   const std::vector<std::string>& want) const {
    const json& v = At(j, key, where);
    if (!v.is_array()) Fail("BBR024", path, where + key + " must be an array of column names");
    std::vector<std::string> out;
    for (const auto& x : v) {
      if (!x.is_string()) break;
      out.push_back(x.get<std::string>());
    }
    if (out.size() != v.size() || out != want)
      Fail("BBR025", path, where + key + " differ from the positional order BBRHFSSData reads");
    return out;
  }
  BBRAxisRange Range(const json& j, const std::string& key, const std::string& where) const {
    const json& r = At(j, key, where);
    const std::string w = where + key + ".";
    BBRAxisRange out;
    out.min = Num(r, "min", w);
    out.max = Num(r, "max", w);
    out.count = Int(r, "count", w);
    if (out.count < 1 || out.max < out.min) Fail("BBR024", path, where + key + " needs min <= max and count >= 1");
    return out;
  }
  void Expect(const std::string& got, const std::string& want, const std::string& field) const {
    if (got != want) Fail("BBR025", path, field + " is \"" + got + "\"; BBRsim implements only \"" + want + "\"");
  }
};

}  // namespace

G4bool BBRCrossSection::Contains(G4double y, G4double z, G4double relTol) const {
  if (shape == kRectangle)
    return std::abs(y) <= yHalf_m * (1. + relTol) && std::abs(z) <= zHalf_m * (1. + relTol);
  return y * y + z * z <= radius_m * radius_m * (1. + relTol);
}

std::vector<std::pair<G4double, G4double>> BBRCrossSection::Boundary(G4double relTol, G4int nDisc) const {
  std::vector<std::pair<G4double, G4double>> out;
  if (shape == kRectangle) {
    const double y = yHalf_m * (1. + relTol), z = zHalf_m * (1. + relTol);
    for (double sy : {-1., 0., 1.})
      for (double sz : {-1., 0., 1.})
        if (sy != 0. || sz != 0.) out.emplace_back(sy * y, sz * z);
  } else {
    const double r = radius_m * std::sqrt(1. + relTol);
    for (G4int i = 0; i < nDisc; ++i) {
      const double a = 2. * kPi * i / nDisc;
      out.emplace_back(r * std::cos(a), r * std::sin(a));
    }
  }
  return out;
}

BBRDatasetSidecar BBRDatasetSidecar::Load(const std::string& waveguidesDir, const std::string& datasetId,
                                          const std::string& dirStem, G4double dirFreqGHz) {
  std::error_code ec;
  const auto abs = std::filesystem::absolute(std::filesystem::path(waveguidesDir) / (dirStem + ".dataset.json"), ec);
  const std::string path = ec ? waveguidesDir + "/" + dirStem + ".dataset.json" : abs.string();
  std::ifstream f(path, std::ios::binary);
  if (!f)
    Fail("BBR024", path, "cannot be opened. Every HFSS dataset needs this sidecar (schema 1.x) beside its two "
                         "_Ephi=N directories; see validation/README.md, Dataset sidecars.");
  std::stringstream ss;
  ss << f.rdbuf();
  return Parse(ss.str(), path, datasetId, dirStem, dirFreqGHz);
}

BBRDatasetSidecar BBRDatasetSidecar::Parse(const std::string& text, const std::string& path,
                                           const std::string& datasetId, const std::string& dirStem,
                                           G4double dirFreqGHz) {
  json j;
  try {
    j = json::parse(text);
  } catch (const json::exception& e) {
    Fail("BBR024", path, std::string("not valid JSON: ") + e.what());
  }
  if (!j.is_object()) Fail("BBR024", path, "the top level must be a JSON object");
  const Reader r{path};
  BBRDatasetSidecar sc;
  sc.path = path;

  // F1-F3: schema, identity, frequency.
  const std::string schema = r.Str(j, "schema_version", "");
  if (schema != "1" && schema.rfind("1.", 0) != 0) Fail("BBR024", path, "schema_version " + schema + ": BBRsim reads 1.x");
  sc.datasetId = r.Str(j, "dataset_id", "");
  if (sc.datasetId != datasetId)
    Fail("BBR024", path, "dataset_id \"" + sc.datasetId + "\", but the directory and crack volume say \"" + datasetId + "\"");
  const std::string prefix = datasetId + "_";
  if (dirStem.compare(0, prefix.size(), prefix) != 0) Fail("BBR024", path, "directory stem " + dirStem + " does not start with " + prefix);
  sc.frequencyLabel = r.Str(j, "frequency_label", "");
  if (sc.frequencyLabel != dirStem.substr(prefix.size()))
    Fail("BBR024", path, "frequency_label \"" + sc.frequencyLabel + "\", but the directory says \"" + dirStem.substr(prefix.size()) + "\"");
  sc.frequencyGHz = r.Num(j, "frequency_ghz", "");
  const double fromLabel = BBRHFSSData::ParseFrequencyGHz(sc.frequencyLabel);
  if (!(dirFreqGHz > 0.) || !RelNear(sc.frequencyGHz, dirFreqGHz, 1e-3) || !RelNear(fromLabel, dirFreqGHz, 1e-3))
    Fail("BBR024", path, "frequency_ghz " + Fmt(sc.frequencyGHz) + " disagrees with the directory frequency " + Fmt(dirFreqGHz) + " GHz (0.1 %)");

  // F4: frames.
  const json& fr = r.At(j, "frames", "");
  const json& ax = r.At(fr, "hfss_global_axes_in_canonical", "frames.");
  const std::string wa = "frames.hfss_global_axes_in_canonical.";
  const Vec hx = r.Vec3(ax, "x", wa), hy = r.Vec3(ax, "y", wa), hz = r.Vec3(ax, "z", wa);
  if (std::abs(Dot(hx, hx) - 1) > kAxisTol || std::abs(Dot(hy, hy) - 1) > kAxisTol || std::abs(Dot(hz, hz) - 1) > kAxisTol ||
      std::abs(Dot(hx, hy)) > kAxisTol || std::abs(Dot(hy, hz)) > kAxisTol || std::abs(Dot(hz, hx)) > kAxisTol ||
      !Near(Cross(hx, hy), hz))
    Fail("BBR025", path, "frames.hfss_global_axes_in_canonical is not a right-handed orthonormal basis");
  const auto canon = [&](const Vec& v) {
    return Vec{v[0] * hx[0] + v[1] * hy[0] + v[2] * hz[0], v[0] * hx[1] + v[1] * hy[1] + v[2] * hz[1],
               v[0] * hx[2] + v[1] * hy[2] + v[2] * hz[2]};
  };
  const json& ecs = r.At(fr, "exit_cs", "frames.");
  const Vec ex = r.Vec3(ecs, "x", "frames.exit_cs."), ey = r.Vec3(ecs, "y", "frames.exit_cs."), ez = r.Vec3(ecs, "z", "frames.exit_cs.");
  if (!Near(Cross(ex, ey), ez)) Fail("BBR025", path, "frames.exit_cs.z " + Fmt(ez) + " is not x cross y " + Fmt(Cross(ex, ey)));
  const json& ecc = r.At(fr, "exit_cs_axes_in_canonical", "frames.");
  const std::string wc = "frames.exit_cs_axes_in_canonical.";
  const Vec cx = r.Vec3(ecc, "x", wc), cy = r.Vec3(ecc, "y", wc), cz = r.Vec3(ecc, "z", wc);
  if (!Near(canon(ex), cx) || !Near(canon(ey), cy) || !Near(canon(ez), cz))
    Fail("BBR025", path, "frames.exit_cs_axes_in_canonical disagrees with exit_cs mapped through hfss_global_axes_in_canonical");
  if (!Near(canon(r.Vec3(fr, "exit_outward_normal_global", "frames.")), kP) ||
      !Near(canon(r.Vec3(fr, "entrance_outward_normal_global", "frames.")), Vec{-1, 0, 0}))
    Fail("BBR025", path, "the exit face's outward normal must map to +p and the entrance face's to -p");
  if (!Near(hx, kHfssX) || !Near(hy, kHfssY) || !Near(hz, kHfssZ) || !Near(cx, kP) || !Near(cy, kL) || !Near(cz, kG))
    Fail("BBR025", path, "the frames differ from the one the sampler implements (HFSS X = -g, Y = +l, Z = +p; "
                         "exit CS = (p, l, g)); another frame needs a per-dataset transform, not implemented yet");

  // F5: excitation.
  const json& exc = r.At(j, "excitation", "");
  r.Expect(r.Str(exc, "coordinate_system", "excitation."), "global", "excitation.coordinate_system");
  r.Expect(r.Str(exc, "incidence_convention", "excitation."), "arrival_direction", "excitation.incidence_convention");
  if (std::abs(r.Num(exc, "normal_entry_theta_deg", "excitation.") - 180.) > 1e-9)
    Fail("BBR025", path, "excitation.normal_entry_theta_deg must be 180 (k = -r: IWaveTheta = 180 deg is normal entry)");
  r.Expect(r.Str(exc, "polarization_convention", "excitation."), kPolarization, "excitation.polarization_convention");
  sc.incidentPhiDeg = r.Nums(exc, "incident_phi_deg", "excitation.");
  sc.incidentThetaDeg = r.Nums(exc, "incident_theta_deg", "excitation.");

  // F6: far field.
  const json& ff = r.At(j, "far_field", "");
  r.Expect(r.Str(ff, "coordinate_system", "far_field."), "exit_cs", "far_field.coordinate_system");
  r.Expect(r.Str(ff, "definition", "far_field."), "Theta-Phi", "far_field.definition");
  r.Expect(r.Str(ff, "component_basis", "far_field."), "spherical_in_exit_cs", "far_field.component_basis");
  sc.farFieldTheta = r.Range(ff, "theta_deg", "far_field.");
  sc.farFieldPhi = r.Range(ff, "phi_deg", "far_field.");
  sc.farFieldPointsPerKey = r.Int(ff, "points_per_key", "far_field.");
  sc.farFieldColumns = r.Columns(ff, "columns", "far_field.", kFarFieldColumns);

  // F7: exit field.
  const json& xf = r.At(j, "exit_field", "");
  r.Expect(r.Str(xf, "coordinate_system", "exit_field."), "exit_cs", "exit_field.coordinate_system");
  if (!r.Bool(xf, "points_in_si", "exit_field.")) Fail("BBR025", path, "exit_field.points_in_si must be true (X, Y, Z in metres)");
  r.Expect(r.Str(xf, "plane", "exit_field."), "x_e=0", "exit_field.plane");
  const bool inRef = r.Bool(xf, "field_in_ref_cs", "exit_field.");
  if (r.Str(xf, "field_components_frame", "exit_field.") != (inRef ? "exit_cs" : "hfss_global"))
    Fail("BBR025", path, "exit_field.field_in_ref_cs and field_components_frame contradict each other");
  const json& grid = r.At(xf, "grid", "exit_field.");
  sc.exitY = r.Range(grid, "y_e", "exit_field.grid.");
  sc.exitZ = r.Range(grid, "z_e", "exit_field.grid.");
  sc.exitPointsPerKey = r.Int(xf, "points_per_key_retained", "exit_field.");
  const std::string outside = r.Str(xf, "outside_points", "exit_field.");
  if (outside != "none" && outside != "omitted" && outside != "zero")
    Fail("BBR024", path, "exit_field.outside_points must be none, omitted or zero");
  sc.outsidePoints = outside;
  sc.exitFieldColumns = r.Columns(xf, "columns", "exit_field.", kExitColumns);

  // F8: cross-section.
  const json& cs = r.At(xf, "cross_section", "exit_field.");
  const std::string wcs = "exit_field.cross_section.";
  const std::string shape = r.Str(cs, "shape", wcs);
  if (shape == "rectangle") {
    sc.crossSection.shape = BBRCrossSection::kRectangle;
    sc.crossSection.yHalf_m = r.Num(cs, "y_e_half_m", wcs);
    sc.crossSection.zHalf_m = r.Num(cs, "z_e_half_m", wcs);
    if (!(sc.crossSection.yHalf_m > 0.) || !(sc.crossSection.zHalf_m > 0.)) Fail("BBR024", path, "rectangle half-extents must be > 0");
  } else if (shape == "disc") {
    sc.crossSection.shape = BBRCrossSection::kDisc;
    sc.crossSection.radius_m = r.Num(cs, "radius_m", wcs);
    if (!(sc.crossSection.radius_m > 0.)) Fail("BBR024", path, "disc radius_m must be > 0");
  } else if (shape == "polygon") {
    Fail("BBR025", path, "cross_section shape polygon is reserved; BBRsim does not support it yet");
  } else {
    Fail("BBR024", path, "unknown cross_section shape \"" + shape + "\"");
  }

  // F9: transmittance definition.
  const json& tr = r.At(j, "transmittance", "");
  r.Expect(r.Str(tr, "definition", "transmittance."), kTransmittance, "transmittance.definition");
  if (r.Bool(tr, "incoming_includes_cos_theta", "transmittance."))
    Fail("BBR025", path, "transmittance.incoming_includes_cos_theta must be false (IngoingPower = area x Ei^2/(2 c mu0))");

  // F10: symmetry.
  const json& sym = r.At(j, "symmetry", "");
  for (const char* k : {"mirror_l", "mirror_g", "end_to_end"})
    if (!r.Bool(sym, k, "symmetry."))
      Fail("BBR025", path, std::string("symmetry.") + k + " is false; the sampler folds by both transverse mirrors "
                                                          "and serves both ends from one table");
  r.Bool(sym, "rotational", "symmetry.");
  const json& boundaries = r.At(j, "boundaries", "");
  if (!boundaries.is_object()) Fail("BBR024", path, "boundaries must be an object");

  // F11: geometry (compared with the placed solid by BBRCrackLibrary; info only).
  const json& geo = r.At(j, "geometry", "");
  const json& ext = r.At(geo, "extent_mm", "geometry.");
  sc.extentP_mm = r.Num(ext, "p", "geometry.extent_mm.");
  sc.extentL_mm = r.Num(ext, "l", "geometry.extent_mm.");
  sc.extentG_mm = r.Num(ext, "g", "geometry.extent_mm.");
  if (!(sc.extentP_mm > 0.) || !(sc.extentL_mm > 0.) || !(sc.extentG_mm > 0.)) Fail("BBR024", path, "geometry.extent_mm must be > 0");

  // F13: the lowest mode and the polarization-filter limit, re-derived from the section.
  const json& md = r.At(j, "modes", "");
  sc.lowestMode = r.Str(md, "mode", "modes.");
  sc.cutoffGHz = r.Num(md, "cutoff_ghz", "modes.");
  sc.polarizationFilterLimitGHz = r.Num(md, "polarization_filter_limit_ghz", "modes.");
  std::string wantMode;
  double wantCut, wantFilter;
  if (sc.crossSection.shape == BBRCrossSection::kRectangle) {
    const double a = 2 * sc.crossSection.yHalf_m, b = 2 * sc.crossSection.zHalf_m;
    if (a >= b) { wantMode = "TE10"; wantCut = kC / (2 * a) / 1e9; } else { wantMode = "TE01"; wantCut = kC / (2 * b) / 1e9; }
    wantFilter = kC / (2 * b) / 1e9;
  } else {
    wantMode = "TE11";
    wantCut = kX11Prime * kC / (2 * kPi * sc.crossSection.radius_m) / 1e9;
    wantFilter = wantCut;
  }
  if (sc.lowestMode != wantMode || !RelNear(sc.cutoffGHz, wantCut, kModeRelTol) ||
      !RelNear(sc.polarizationFilterLimitGHz, wantFilter, kModeRelTol))
    Fail("BBR025", path, "modes (" + sc.lowestMode + ", " + Fmt(sc.cutoffGHz) + ", " + Fmt(sc.polarizationFilterLimitGHz) +
                         " GHz) disagree with the cross-section (" + wantMode + ", " + Fmt(wantCut) + ", " + Fmt(wantFilter) + " GHz)");

  json inv;
  inv["frames"] = {{"hfss_global_axes_in_canonical", ax}, {"exit_cs_axes_in_canonical", ecc}};
  inv["symmetry"] = sym;
  inv["boundaries"] = boundaries;
  json g = {{"extent_mm", ext}};
  if (geo.contains("shape")) g["shape"] = geo.at("shape");
  inv["geometry"] = g;
  json m = md;
  m.erase("propagating_count");
  m.erase("basis");
  inv["modes"] = m;
  inv["cross_section"] = cs;
  sc.invariant = inv.dump();
  return sc;
}

namespace {
// Objects over the keys both carry, arrays element by element, numbers to 1e-9
// relative, everything else exactly (the frequency-independent comparison).
bool Close(const json& a, const json& b) {
  if (a.is_number() && b.is_number()) {
    const double x = a.get<double>(), y = b.get<double>();
    return std::abs(x - y) <= 1e-9 * std::max(std::abs(x), std::abs(y));
  }
  if (a.is_object() && b.is_object()) {
    for (auto it = a.begin(); it != a.end(); ++it)
      if (b.contains(it.key()) && !Close(it.value(), b.at(it.key()))) return false;
    return true;
  }
  if (a.is_array() && b.is_array()) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (!Close(a[i], b[i])) return false;
    return true;
  }
  return a == b;
}
}  // namespace

G4bool BBRDatasetSidecar::SameInvariant(const BBRDatasetSidecar& other) const {
  if (invariant.empty() || other.invariant.empty()) return false;   // not produced by Parse
  return Close(json::parse(invariant), json::parse(other.invariant));
}

void BBRDatasetSidecar::CheckFitsSolid(const G4VSolid& solid, const std::string& volumeName) const {
  if (solid.Inside(G4ThreeVector()) != kInside)
    Fail("BBR025", path, "the local origin of crack volume " + volumeName + " is not inside its solid; the wrapper "
                         "finds the exit face from the origin along local x (BBSimOpBoundaryProcess.cc)");
  const double inset = 10. * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  int bad = 0;
  G4ThreeVector worst;
  for (const double dir : {1., -1.}) {
    const double half = solid.DistanceToOut(G4ThreeVector(), G4ThreeVector(dir, 0., 0.));
    for (const auto& [y, z] : crossSection.Boundary())
      for (const double sy : {1., -1.})
        for (const double sz : {1., -1.}) {
          const G4ThreeVector p(dir * (half - inset), sy * y * CLHEP::m, sz * z * CLHEP::m);
          if (solid.Inside(p) != kInside) { ++bad; worst = p; }
        }
  }
  if (bad) {
    std::ostringstream ed;
    ed << "the declared exit cross-section does not fit strictly inside crack volume " << volumeName << ": "
       << bad << " boundary sample(s) are not inside, e.g. local " << worst / mm << " mm. Make the Geant4 section "
       << "larger than the HFSS one by a recorded margin (as the 52/102 um test-world gaps are), or check the "
       << "solid's axis (local x must be the propagation axis; a bare G4Tubs has it on z).";
    Fail("BBR025", path, ed.str());
  }
}
