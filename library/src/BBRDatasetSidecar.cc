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
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using json = nlohmann::json;

namespace {

using Vec = std::array<double, 3>;

constexpr double kC = 299792458.0;           // m/s, exact
constexpr double kX11Prime = 1.8411837813;   // first zero of J1': TE11 of a round guide
constexpr double kX01 = 2.4048255577;        // first zero of J0: TM01 of a round guide
constexpr double kX21Prime = 3.0542369282;   // first zero of J2': TE21 of a round guide
constexpr double kListLimitGHz = 20000.;     // the regime list runs to 20 THz, inclusive
constexpr double kPi = 3.14159265358979323846;
constexpr double kAxisTol = 1e-9;
constexpr double kModeRelTol = 1e-6;

const char* const kPolarization = "Ephi=0: E_theta=1; Ephi=1: E_phi=1";
const char* const kTransmittance = "outgoing_power_w / incoming_power_w, first row per key";
const std::vector<std::string> kFarFieldColumns = {"Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta",
                                                   "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"};
// The 2.0 list: headers in the stored (interleaved) order, InterceptFactor last. Read by position:
// Ex = (9, 10), Ey = (11, 12), Ez = (13, 14), InterceptFactor = 15.
const std::vector<std::string> kExitColumns2 = {"Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower",
                                                "IngoingPower", "X", "Y", "Z", "Ex_real", "Ex_imag", "Ey_real",
                                                "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor"};
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

// ---- schema 2.x --------------------------------------------------------------------------------
constexpr double kKeyTol = 1e-9;       // per-key records match the declared keys (degrees)
constexpr double kDivisorRelTol = 1e-6;   // c_k, cell_90 against the closed form
constexpr double kStepRelTol = 1e-6;
constexpr double kCosFloor = 0.01;

std::string KeyStr(double phi, double theta) { return "(phi " + Fmt(phi) + ", theta " + Fmt(theta) + ")"; }

// Delta = (smallest declared theta above 90 deg) - 90, or -1 without a 90 deg entry or anything above it.
double Delta(const std::vector<double>& thetaDeg) {
  bool has90 = false;
  double above = 1e300;
  for (const double t : thetaDeg) {
    if (std::abs(t - 90.) <= kKeyTol) has90 = true;
    else if (t > 90.) above = std::min(above, t);
  }
  return (has90 && above < 1e300) ? above - 90. : -1.;
}

// The records of one per-key list, indexed phi-major then theta, each matched to a declared key
// within kKeyTol; a duplicate, an unmatched record or a gap is BBR024.
std::vector<const json*> MatchKeys(const Reader& r, const json& parent, const std::string& key, const std::string& where,
                                   const std::vector<double>& phis, const std::vector<double>& thetas) {
  const json& list = r.At(parent, key, where);
  const std::string w = where + key;
  if (!list.is_array()) Fail("BBR024", r.path, w + " must be a list of per-key records");
  std::vector<const json*> out(phis.size() * thetas.size(), nullptr);
  for (std::size_t i = 0; i < list.size(); ++i) {
    const std::string wi = w + "[" + std::to_string(i) + "].";
    if (!list[i].is_object()) Fail("BBR024", r.path, w + "[" + std::to_string(i) + "] must be an object");
    const double phi = r.Num(list[i], "phi_deg", wi), theta = r.Num(list[i], "theta_deg", wi);
    std::size_t ip = phis.size(), it = thetas.size();
    for (std::size_t a = 0; a < phis.size(); ++a) if (std::abs(phi - phis[a]) <= kKeyTol) { ip = a; break; }
    for (std::size_t b = 0; b < thetas.size(); ++b) if (std::abs(theta - thetas[b]) <= kKeyTol) { it = b; break; }
    if (ip == phis.size() || it == thetas.size())
      Fail("BBR024", r.path, wi + "key " + KeyStr(phi, theta) + " does not match a declared key of excitation.incident_phi_deg x incident_theta_deg");
    const std::size_t slot = ip * thetas.size() + it;
    if (out[slot]) Fail("BBR024", r.path, w + ": duplicate record for key " + KeyStr(phis[ip], thetas[it]));
    out[slot] = &list[i];
  }
  for (std::size_t ip = 0; ip < phis.size(); ++ip)
    for (std::size_t it = 0; it < thetas.size(); ++it)
      if (!out[ip * thetas.size() + it])
        Fail("BBR024", r.path, w + " has no record for declared key " + KeyStr(phis[ip], thetas[it]));
  return out;
}

// The two polarizations of a per-key record: exactly "0" and "1", as in the field list.
std::array<const json*, 2> Polarizations(const Reader& r, const json& rec, const std::string& where) {
  const json& pols = r.At(rec, "polarizations", where);
  if (!pols.is_object() || pols.size() != 2) Fail("BBR024", r.path, where + "polarizations must be an object with exactly \"0\" and \"1\"");
  return {&r.At(pols, "0", where + "polarizations."), &r.At(pols, "1", where + "polarizations.")};
}

void ParseSchema2(const Reader& r, const json& j, BBRDatasetSidecar& sc) {
  const std::string& path = r.path;
  r.Expect(r.Str(r.At(j, "model", ""), "kind", "model."), "screen_free_abc", "model.kind");
  if (!r.At(j, "conventions", "").is_object()) Fail("BBR024", path, "conventions must be an object");

  // transmittance (OR-034, HF-049): T = OutgoingPower / (IngoingPower * InterceptFactor).
  const json& tr = r.At(j, "transmittance", "");
  const double want90 = BBRDatasetSidecar::Cell90(sc.incidentThetaDeg);
  if (want90 < 0.)
    Fail("BBR025", path, "excitation.incident_theta_deg has no theta = 90 deg entry with a declared theta above it, so "
                         "transmittance.cell_90 cannot be derived from the declared grid (HF-049); BBRsim refuses it");
  const double delta = Delta(sc.incidentThetaDeg);
  r.Expect(r.Str(tr, "definition", "transmittance."), "OutgoingPower / (IngoingPower * InterceptFactor)", "transmittance.definition");
  r.Expect(r.Str(tr, "numerator", "transmittance."), "exit_flux_re_s_dot_n", "transmittance.numerator");
  r.Expect(r.Str(tr, "incoming", "transmittance."), "area_times_intensity_undivided", "transmittance.incoming");
  r.Expect(r.Str(tr, "above_one", "transmittance."), "per_photon_cap", "transmittance.above_one");
  sc.cell90 = r.Num(tr, "cell_90", "transmittance.");
  if (!RelNear(sc.cell90, want90, kDivisorRelTol))
    Fail("BBR025", path, "transmittance.cell_90 is " + Fmt(sc.cell90) + ", but the declared grid's outermost band gives " + Fmt(want90) +
                         " (1e-6 relative; HF-049)");
  const json& as = r.At(tr, "cell_90_assumes", "transmittance.");
  const std::string wa = "transmittance.cell_90_assumes.";
  r.Expect(r.Str(as, "emitter", wa), "uniform_theta", wa + "emitter");   // BBRsim's one emitter (ThermalSurface.cc:133)
  const double floorV = r.Num(as, "cos_floor", wa);
  if (!RelNear(floorV, kCosFloor, 1e-9))
    Fail("BBR025", path, wa + "cos_floor is " + Fmt(floorV) + "; BBRsim implements only " + Fmt(kCosFloor));
  const std::vector<double> band = r.Nums(as, "band_deg", wa);
  if (band.size() != 2) Fail("BBR024", path, wa + "band_deg must be [lo, 90]");
  if (std::abs(band[0] - (90. - delta / 2)) > kKeyTol || std::abs(band[1] - 90.) > kKeyTol)
    Fail("BBR025", path, wa + "band_deg is [" + Fmt(band[0]) + ", " + Fmt(band[1]) + "]; the declared grid gives [" +
                         Fmt(90. - delta / 2) + ", 90]");

  const auto& phis = sc.incidentPhiDeg;
  const auto& thetas = sc.incidentThetaDeg;
  const auto cKeyRecs = MatchKeys(r, tr, "c_key", "transmittance.", phis, thetas);
  const auto perKeyRecs = MatchKeys(r, tr, "per_key", "transmittance.", phis, thetas);
  const json& xf = j.at("exit_field");
  const auto scaleRecs = MatchKeys(r, r.At(xf, "scale", "exit_field."), "per_key", "exit_field.scale.", phis, thetas);
  const auto ffRecs = MatchKeys(r, r.At(j, "far_field", ""), "ff_over_flux", "far_field.", phis, thetas);

  sc.keys.clear();
  for (std::size_t ip = 0; ip < phis.size(); ++ip)
    for (std::size_t it = 0; it < thetas.size(); ++it) {
      const std::size_t k = ip * thetas.size() + it;
      const std::string ks = KeyStr(phis[ip], thetas[it]);
      BBRKeyRecord rec;
      rec.phiDeg = phis[ip];
      rec.thetaDeg = thetas[it];
      rec.cK = r.Num(*cKeyRecs[k], "c_k", "transmittance.c_key " + ks + ".");
      const double wantC = BBRDatasetSidecar::InterceptFactor(rec.thetaDeg, thetas);
      if (!RelNear(rec.cK, wantC, kDivisorRelTol))
        Fail("BBR025", path, "transmittance.c_key " + ks + " c_k is " + Fmt(rec.cK) + ", but the declared grid gives " + Fmt(wantC) +
                             " (|cos theta|, or cell_90 at 90 deg; 1e-6 relative)");
      const json& pk = *perKeyRecs[k];
      const std::string wp = "transmittance.per_key " + ks + ".";
      const double pkC = r.Num(pk, "c_k", wp);
      if (!RelNear(pkC, rec.cK, kDivisorRelTol))
        Fail("BBR024", path, wp + "c_k " + Fmt(pkC) + " disagrees with transmittance.c_key " + Fmt(rec.cK));
      const json& sk = *scaleRecs[k];
      const std::string ws = "exit_field.scale.per_key " + ks + ".";
      rec.w0 = r.Num(sk, "w0", ws);
      rec.w1 = r.Num(sk, "w1", ws);
      if (!(rec.w0 > 0.) || !(rec.w1 > 0.)) Fail("BBR024", path, ws + "w0 and w1 must be > 0 (w0 " + Fmt(rec.w0) + ", w1 " + Fmt(rec.w1) + ")");
      rec.dark0 = r.Bool(sk, "dark0", ws);
      rec.dark1 = r.Bool(sk, "dark1", ws);
      const auto pols = Polarizations(r, pk, wp);
      const bool dark[2] = {rec.dark0, rec.dark1};
      double* tf[2] = {&rec.tFlux0, &rec.tFlux1};
      for (int p = 0; p < 2; ++p) {
        const std::string wpp = wp + "polarizations." + std::to_string(p) + ".";
        *tf[p] = r.Num(*pols[p], "t_flux", wpp);
        r.Num(*pols[p], "outgoing_power_abs_re_s_w", wpp);
        // OR-038: a dark basis has t_flux exactly 0; the sidecar must not contradict itself.
        if (dark[p] && *tf[p] != 0.)
          Fail("BBR024", path, "key " + ks + " polarization " + std::to_string(p) + " is dark (exit_field.scale.per_key dark" +
                               std::to_string(p) + "), but transmittance.per_key t_flux is " + Fmt(*tf[p]) +
                               "; a dark basis needs t_flux exactly 0 (OR-038)");
      }
      const auto ffp = Polarizations(r, *ffRecs[k], "far_field.ff_over_flux " + ks + ".");
      for (int p = 0; p < 2; ++p) {
        const json& v = *ffp[p];
        const std::string wf = "far_field.ff_over_flux " + ks + " polarizations." + std::to_string(p);
        if (v.is_null()) {
          if (!dark[p]) Fail("BBR024", path, wf + " is null, but polarization " + std::to_string(p) + " of key " + ks + " is not dark");
        } else if (!v.is_number() || !std::isfinite(v.get<double>())) {
          Fail("BBR024", path, wf + " must be a finite number, or null for a dark basis");
        }
      }
      sc.keys.push_back(rec);
    }

  // exit field: the sampler, quadrature and normalization the reader implements, and the fitted lattice.
  const bool rectangle = sc.crossSection.shape == BBRCrossSection::kRectangle;
  r.Expect(r.Str(xf, "quadrature", "exit_field."), rectangle ? "trapezoid" : "uniform_disc", "exit_field.quadrature");
  r.Expect(r.Str(xf, "sampler", "exit_field."), "abs_e_t_squared_times_quadrature", "exit_field.sampler");
  const json& comps = r.At(xf, "tangential_components", "exit_field.");
  if (!comps.is_array()) Fail("BBR024", path, "exit_field.tangential_components must be an array");
  if (!(comps.size() == 2 && comps[0] == "Ex" && comps[1] == "Ey"))
    Fail("BBR025", path, "exit_field.tangential_components must be [\"Ex\", \"Ey\"]; the sampler reads the exit-face tangentials only");
  r.Expect(r.Str(r.At(xf, "scale", "exit_field."), "rule", "exit_field.scale."),
           "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w", "exit_field.scale.rule");
  const json& co = r.At(xf, "column_order_check", "exit_field.");
  if (!co.is_object()) Fail("BBR024", path, "exit_field.column_order_check must be an object");
  r.Expect(r.Str(co, "rule", "exit_field.column_order_check."), "reference_exporter_interleaved_v1", "exit_field.column_order_check.rule");
  r.Str(co, "method", "exit_field.column_order_check.");
  const json& ref = r.At(co, "reference", "exit_field.column_order_check.");
  if (!ref.is_object()) Fail("BBR024", path, "exit_field.column_order_check.reference must be an object");
  if (!r.Bool(ref, "passed", "exit_field.column_order_check.reference."))
    Fail("BBR025", path, "exit_field.column_order_check.reference.passed is false: the exporter's column order was not validated");
  const double planeTol = r.Num(xf, "plane_tolerance_m", "exit_field.");
  if (!RelNear(planeTol, 1e-12, 1e-9))
    Fail("BBR025", path, "exit_field.plane_tolerance_m is " + Fmt(planeTol) + "; BBRsim implements only 1e-12");
  const json& grid = xf.at("grid");
  for (const char* axis : {"y_e", "z_e"}) {
    BBRAxisRange& a = std::string(axis) == "y_e" ? sc.exitY : sc.exitZ;
    const std::string wg = std::string("exit_field.grid.") + axis;
    const json& ga = grid.at(axis);
    if (rectangle && a.count < 2)
      Fail("BBR024", path, wg + " needs count >= 2 on a rectangle: neither the lattice step nor the trapezoid weights are defined on one node");
    a.step = r.Num(ga, "step", wg + ".");
    if (!(a.step > 0.)) Fail("BBR024", path, wg + ".step must be > 0");
    if (a.count >= 2) {
      const double want = (a.max - a.min) / (a.count - 1);
      if (!RelNear(a.step, want, kStepRelTol))
        Fail("BBR024", path, wg + ".step is " + Fmt(a.step) + ", but (max - min)/(count - 1) is " + Fmt(want));
    }
  }

  // cross term and far field.
  const json& ct = r.At(j, "cross_term", "");
  r.Expect(r.Str(ct, "density", "cross_term."), "e_t_quadrature", "cross_term.density");
  r.Str(ct, "definition", "cross_term.");
  sc.rhoValidBelowGHz = r.Num(ct, "valid_below_ghz", "cross_term.");
  r.Expect(r.Str(ct, "above_limit", "cross_term."), "warn_once", "cross_term.above_limit");
  r.Expect(r.Str(r.At(j, "far_field", ""), "normalization", "far_field."), "pattern_only", "far_field.normalization");
}

}  // namespace

std::vector<std::pair<std::string, G4double>> BBRDatasetSidecar::RegimeCutoffs(const BBRCrossSection& s) {
  std::vector<std::pair<std::string, G4double>> out;
  if (s.shape == BBRCrossSection::kDisc) {
    // Closed-PEC round guide: f = x c / (2 pi R) with the zeros x'11 (TE11), x01 (TM01) and x'21 (TE21).
    const double k = kC / (2 * kPi * s.radius_m) / 1e9;
    out.emplace_back("TE11", kX11Prime * k);
    out.emplace_back("TM01", kX01 * k);
    out.emplace_back("TE21", kX21Prime * k);
    return out;
  }
  // Closed-PEC rectangle, a = 2 y_e_half along l, b = 2 z_e_half along g: TE10 at c/(2a) and the
  // gap family TE0n at n c/(2b) up to the list limit, inclusive.
  const double a = 2 * s.yHalf_m, b = 2 * s.zHalf_m;
  out.emplace_back("TE10", kC / (2 * a) / 1e9);
  for (int n = 1; n * kC / (2 * b) / 1e9 <= kListLimitGHz + 1e-9; ++n)
    out.emplace_back(n < 10 ? "TE0" + std::to_string(n) : "TE0," + std::to_string(n), n * kC / (2 * b) / 1e9);
  std::stable_sort(out.begin(), out.end(), [](const auto& x, const auto& y) { return x.second < y.second; });
  return out;
}

G4double BBRDatasetSidecar::Cell90(const std::vector<G4double>& thetaDeg) {
  const double delta = Delta(thetaDeg);
  if (delta < 0.) return -1.;
  const double a0 = kPi / 2 - delta / 2 * kPi / 180., af = std::acos(kCosFloor), L = kPi / 2 - a0;
  const auto prim = [](double a) { return std::log(1. / std::cos(a) + std::tan(a)); };   // integral of sec
  const double I = af > a0 ? (prim(af) - prim(a0)) + (kPi / 2 - af) / kCosFloor : (kPi / 2 - a0) / kCosFloor;
  return L / I;
}

G4double BBRDatasetSidecar::InterceptFactor(G4double thetaDeg, const std::vector<G4double>& thetaList) {
  if (std::abs(thetaDeg - 90.) <= kKeyTol) return Cell90(thetaList);
  return std::abs(std::cos(thetaDeg * kPi / 180.));
}

const BBRKeyRecord* BBRDatasetSidecar::Key(G4double phiDeg, G4double thetaDeg) const {
  for (const auto& k : keys)
    if (std::abs(k.phiDeg - phiDeg) <= kKeyTol && std::abs(k.thetaDeg - thetaDeg) <= kKeyTol) return &k;
  return nullptr;
}

G4bool BBRCrossSection::Contains(G4double y, G4double z, G4double relTol) const {
  if (shape == kRectangle)
    return std::abs(y) <= yHalf_m * (1. + relTol) && std::abs(z) <= zHalf_m * (1. + relTol);
  return y * y + z * z <= radius_m * radius_m * (1. + relTol);
}

std::vector<std::pair<G4double, G4double>> BBRCrossSection::Boundary(G4double relTol, G4int nDisc,
                                                                     G4double offset_m) const {
  std::vector<std::pair<G4double, G4double>> out;
  if (shape == kRectangle) {
    const double y = yHalf_m * (1. + relTol) + offset_m, z = zHalf_m * (1. + relTol) + offset_m;
    for (double sy : {-1., 0., 1.})
      for (double sz : {-1., 0., 1.})
        if (sy != 0. || sz != 0.) out.emplace_back(sy * y, sz * z);
  } else {
    const double r = radius_m * std::sqrt(1. + relTol) + offset_m;
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
  const bool v1 = schema == "1" || schema.rfind("1.", 0) == 0, v2 = schema == "2" || schema.rfind("2.", 0) == 0;
  if (!v1 && !v2) Fail("BBR024", path, "schema_version " + schema + ": BBRsim reads 1.x and 2.x");
  sc.schemaMajor = v2 ? 2 : 1;
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
  sc.exitFieldColumns = r.Columns(xf, "columns", "exit_field.", sc.schemaMajor == 2 ? kExitColumns2 : kExitColumns);

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

  // F9: transmittance definition. Under 2.x the transmittance, exit-field, cross-term and per-key
  // blocks replace it (OR-034, HF-040, HF-042, HF-047).
  if (sc.schemaMajor == 1) {
    const json& tr = r.At(j, "transmittance", "");
    r.Expect(r.Str(tr, "definition", "transmittance."), kTransmittance, "transmittance.definition");
    if (r.Bool(tr, "incoming_includes_cos_theta", "transmittance."))
      Fail("BBR025", path, "transmittance.incoming_includes_cos_theta must be false (IngoingPower = area x Ei^2/(2 c mu0))");
  } else {
    ParseSchema2(r, j, sc);
  }

  // F10: symmetry.
  const json& sym = r.At(j, "symmetry", "");
  for (const char* k : {"mirror_l", "mirror_g", "end_to_end"})
    if (!r.Bool(sym, k, "symmetry."))
      Fail("BBR025", path, std::string("symmetry.") + k + " is false; the sampler folds by both transverse mirrors "
                                                          "and serves both ends from one table");
  r.Bool(sym, "rotational", "symmetry.");
  // The boundary model the transport implements: radiation end faces and PEC walls. A
  // non-transmission is booked as specular reflection with no wall loss, and the mode formulas
  // assume a closed PEC guide (CDX-003).
  const json& boundaries = r.At(j, "boundaries", "");
  if (!boundaries.is_object()) Fail("BBR024", path, "boundaries must be an object");
  r.Expect(r.Str(boundaries, "entrance", "boundaries."), "radiation", "boundaries.entrance");
  r.Expect(r.Str(boundaries, "exit", "boundaries."), "radiation", "boundaries.exit");
  r.Expect(r.Str(boundaries, "walls", "boundaries."), "PEC", "boundaries.walls");

  // The gap must be vacuum (CDX-029): provenance.material, the producer's AEDT material name, must read
  // "vacuum" in any case when it is a string. Under 1.x it may be null or absent (the legacy sidecars
  // carry null); under 2.x it is required (OR-039), so null, an absent key or no provenance is BBR024.
  const bool hasMaterial = j.contains("provenance") && j.at("provenance").is_object() && j.at("provenance").contains("material");
  if (sc.schemaMajor == 2 && (!hasMaterial || j.at("provenance").at("material").is_null()))
    Fail("BBR024", path, "provenance.material is missing or null; schema 2.x requires it (the gap must be vacuum)");
  if (hasMaterial) {
    const json& mat = j.at("provenance").at("material");
    if (mat.is_string()) {
      std::string lower = mat.get<std::string>();
      for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      if (lower != "vacuum")
        Fail("BBR025", path, "provenance.material is \"" + mat.get<std::string>() + "\"; the gap must be vacuum");
    } else if (!mat.is_null()) {
      Fail("BBR024", path, "provenance.material must be a string or null");
    }
  }

  // F11: geometry. p is compared with the placed solid's length in F12 (CheckFitsSolid); for a rectangle
  // l and g must equal the full widths of the declared exit section (CDX-028). A disc's l and g are not
  // compared with its diameter here: the round gap's extent_mm is the HFSS bounding box.
  const json& geo = r.At(j, "geometry", "");
  const json& ext = r.At(geo, "extent_mm", "geometry.");
  sc.extentP_mm = r.Num(ext, "p", "geometry.extent_mm.");
  sc.extentL_mm = r.Num(ext, "l", "geometry.extent_mm.");
  sc.extentG_mm = r.Num(ext, "g", "geometry.extent_mm.");
  if (!(sc.extentP_mm > 0.) || !(sc.extentL_mm > 0.) || !(sc.extentG_mm > 0.)) Fail("BBR024", path, "geometry.extent_mm must be > 0");
  if (sc.crossSection.shape == BBRCrossSection::kRectangle) {
    const struct { const char* name; double extent_mm; const char* half; double width_mm; } axes[] = {
      {"l", sc.extentL_mm, "y_e_half_m", 2e3 * sc.crossSection.yHalf_m},
      {"g", sc.extentG_mm, "z_e_half_m", 2e3 * sc.crossSection.zHalf_m}};
    for (const auto& a : axes)
      if (std::abs(a.extent_mm - a.width_mm) > 1e-6 * a.width_mm) {
        std::ostringstream ed;
        ed << std::setprecision(17) << "geometry.extent_mm." << a.name << " = " << a.extent_mm << " mm differs from the full width "
           << a.width_mm << " mm of exit_field.cross_section (2 x " << a.half << ") by more than 1e-6 relative";
        Fail("BBR024", path, ed.str());
      }
  }

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

  // OR-031: the cutoffs a photon can cross to reach a table. Under 2.x the declared list, re-derived
  // entry by entry; under 1.x only the lowest mode, so the schema-1 guard is unchanged.
  if (sc.schemaMajor == 1) {
    sc.regimeCutoffs = {{sc.lowestMode, sc.cutoffGHz}};
  } else {
    const json& list = r.At(md, "regime_cutoffs_ghz", "modes.");
    if (!list.is_array()) Fail("BBR024", path, "modes.regime_cutoffs_ghz must be an array of {mode, cutoff_ghz}");
    const auto want = RegimeCutoffs(sc.crossSection);
    for (const auto& e : list) {
      if (!e.is_object()) Fail("BBR024", path, "modes.regime_cutoffs_ghz entries must be objects {mode, cutoff_ghz}");
      sc.regimeCutoffs.emplace_back(r.Str(e, "mode", "modes.regime_cutoffs_ghz[]."), r.Num(e, "cutoff_ghz", "modes.regime_cutoffs_ghz[]."));
    }
    const std::size_t n = std::max(sc.regimeCutoffs.size(), want.size());
    for (std::size_t i = 0; i < n; ++i) {
      const std::string at = "modes.regime_cutoffs_ghz[" + std::to_string(i) + "]";
      if (i >= sc.regimeCutoffs.size())
        Fail("BBR025", path, at + " is missing: the section gives " + want[i].first + " at " + Fmt(want[i].second) + " GHz (" +
                             std::to_string(want.size()) + " entries, the sidecar declares " + std::to_string(sc.regimeCutoffs.size()) + ")");
      if (i >= want.size())
        Fail("BBR025", path, at + " is " + sc.regimeCutoffs[i].first + " at " + Fmt(sc.regimeCutoffs[i].second) + " GHz, but the section gives only " +
                             std::to_string(want.size()) + " entries");
      if (sc.regimeCutoffs[i].first != want[i].first || !RelNear(sc.regimeCutoffs[i].second, want[i].second, kModeRelTol))
        Fail("BBR025", path, at + " is " + sc.regimeCutoffs[i].first + " at " + Fmt(sc.regimeCutoffs[i].second) + " GHz, but the section gives " +
                             want[i].first + " at " + Fmt(want[i].second) + " GHz (1e-6 relative)");
    }
    if (!RelNear(sc.rhoValidBelowGHz, sc.polarizationFilterLimitGHz, kModeRelTol))
      Fail("BBR025", path, "cross_term.valid_below_ghz is " + Fmt(sc.rhoValidBelowGHz) + ", but modes.polarization_filter_limit_ghz is " +
                           Fmt(sc.polarizationFilterLimitGHz) + " GHz; the cross term is valid below that limit (1e-6 relative)");
  }

  json inv;
  inv["schema"] = sc.schemaMajor;
  inv["excitation"] = {{"incident_phi_deg", exc.at("incident_phi_deg")}, {"incident_theta_deg", exc.at("incident_theta_deg")}};
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

  // Recorded only, never checked: an absent field or parent object is null.
  const auto field = [&j](std::initializer_list<const char*> keys) -> json {
    const json* v = &j;
    for (const char* k : keys) {
      if (!v->is_object() || !v->contains(k)) return nullptr;
      v = &v->at(k);
    }
    return *v;
  };
  json rec;
  rec["provenance"] = field({"provenance"});
  rec["boundaries"] = field({"boundaries"});
  rec["geometry"] = {{"shape", field({"geometry", "shape"})}, {"bounding_box_mm", field({"geometry", "bounding_box_mm"})}};
  rec["symmetry"] = {{"rotational", field({"symmetry", "rotational"})}};
  rec["frames"] = {{"pose_rule", field({"frames", "pose_rule"})},
                   {"exit_cs_origin_mm_global", field({"frames", "exit_cs", "origin_mm_global"})}};
  rec["excitation"] = {{"origin_mm_global", field({"excitation", "origin_mm_global"})}};
  sc.recorded = rec.dump(-1, ' ', false, json::error_handler_t::replace);
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

std::vector<std::string> BBRDatasetSidecar::InvariantDiff(const BBRDatasetSidecar& other) const {
  std::vector<std::string> differ;
  if (invariant.empty() || other.invariant.empty()) return differ;   // not produced by Parse
  const json a = json::parse(invariant), b = json::parse(other.invariant);
  for (const char* k : {"schema", "excitation", "frames", "symmetry", "boundaries", "geometry", "modes", "cross_section"})
    if (a.contains(k) && b.contains(k) && !Close(a.at(k), b.at(k))) differ.emplace_back(k);
  return differ;
}

G4bool BBRDatasetSidecar::SameInvariant(const BBRDatasetSidecar& other) const {
  if (invariant.empty() || other.invariant.empty()) return false;   // not produced by Parse
  return InvariantDiff(other).empty();
}

void BBRDatasetSidecar::CheckFitsSolid(const G4VSolid& solid, const std::string& volumeName) const {
  if (solid.Inside(G4ThreeVector()) != kInside)
    Fail("BBR025", path, "the local origin of crack volume " + volumeName + " is not inside its solid; the wrapper "
                         "finds the exit face from the origin along local x (BBSimOpBoundaryProcess.cc)");
  // F12 (OR-020): the solid must equal the declared section. On the sampling plane just inside each
  // end face, boundary points kClampInset_m inside the declared section must be inside the solid (the
  // solid is at most that much smaller), and points kClampInset_m beyond the section inflated by the
  // 1e-6 relative containment tolerance must not be (the solid is at most that much larger). The
  // clamp at load moves exit points to kClampInset_m inside the section, so both bounds are what a
  // clamped point needs.
  const double inset = 20. * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  const auto inward = crossSection.Boundary(0., 360, -BBRCrossSection::kClampInset_m);
  const auto outward = crossSection.Boundary(1e-6, 360, BBRCrossSection::kClampInset_m);
  int small = 0, large = 0;
  G4ThreeVector worst;
  for (const double dir : {1., -1.}) {
    const double half = solid.DistanceToOut(G4ThreeVector(), G4ThreeVector(dir, 0., 0.));
    for (const double sy : {1., -1.})
      for (const double sz : {1., -1.}) {
        for (const auto& [y, z] : inward) {
          const G4ThreeVector p(dir * (half - inset), sy * y * CLHEP::m, sz * z * CLHEP::m);
          if (solid.Inside(p) != kInside) { ++small; worst = p; }
        }
        for (const auto& [y, z] : outward) {
          const G4ThreeVector p(dir * (half - inset), sy * y * CLHEP::m, sz * z * CLHEP::m);
          if (solid.Inside(p) == kInside) { ++large; worst = p; }
        }
      }
  }
  if (small || large) {
    std::ostringstream ed;
    ed << "the declared exit cross-section does not equal crack volume " << volumeName << "'s section: "
       << small << " sample(s) just inside the declared boundary are not inside the solid and " << large
       << " just outside it are inside, e.g. local " << worst / mm << " mm. The solid is "
       << (small && large ? "smaller than the HFSS section in places and larger than it in others"
                          : small ? "smaller than the HFSS section" : "larger than the HFSS section")
       << "; it must equal it (OR-020) within 1 nm, or 1 nm plus 1e-6 relative where larger, with local x "
       << "the propagation axis (a bare G4Tubs has it on z).";
    Fail("BBR025", path, ed.str());
  }
  // Length (CDX-028, OR-020): the solid's extent along local x from the origin equals extent_mm.p within 1 nm,
  // or 1 nm plus 1e-6 relative where longer.
  const double length_mm = (solid.DistanceToOut(G4ThreeVector(), G4ThreeVector(1., 0., 0.)) +
                            solid.DistanceToOut(G4ThreeVector(), G4ThreeVector(-1., 0., 0.))) / mm;
  const double nm_mm = 1e-6;
  if (length_mm < extentP_mm - nm_mm || length_mm > extentP_mm * (1. + 1e-6) + nm_mm) {
    std::ostringstream ed;
    ed << std::setprecision(17) << "the crack volume's length " << length_mm << " mm differs from the HFSS model's "
       << extentP_mm << " mm (OR-020) for crack volume " << volumeName << ": the solid's length along local x must "
       << "equal geometry.extent_mm.p within 1 nm, or 1 nm plus 1e-6 relative where longer";
    Fail("BBR025", path, ed.str());
  }
}
