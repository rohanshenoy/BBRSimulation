#include "BBRHFSSData.hh"
#include "BBRDatasetSidecar.hh"

#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> SplitCSV(std::string line)
{
  // The shipped CSVs have CRLF line ends and std::getline keeps the '\r';
  // drop it so the last field parses in full (Num rejects trailing bytes).
  if (!line.empty() && line.back() == '\r') line.pop_back();
  std::vector<std::string> v;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, ',')) v.push_back(tok);
  return v;
}

// Stable floating-point map key: round to nearest 0.01 degree.
G4double RoundDeg(G4double d) { return std::round(d * 100.) / 100.; }

// Far-field grid directions whose component along the exit-face normal is
// below this carry no flux through the face and are never sampled. The HFSS
// Phi sweep runs -90.0 .. +89.67 deg, so Phi = -90 (a direction lying IN the
// face plane, cos = 6e-17) is on the grid; sampling it left a transmitted
// photon skating along the exit face into the crack wall.
constexpr G4double kMinNormalComponent = 1e-6;

std::pair<G4double, G4double> MakeKey(G4double phi, G4double theta)
{
  return {RoundDeg(phi), RoundDeg(theta)};
}

// std::stod that turns a bad or non-finite field into BBR013 instead of a
// C++ exception (which would terminate a worker thread).
G4double Num(const std::string& s, const G4String& path, const char* what)
{
  try {
    std::size_t used = 0;
    const G4double v = std::stod(s, &used);
    if (used != s.size() || !std::isfinite(v)) throw std::invalid_argument(s);
    return v;
  } catch (const std::exception&) {
    G4Exception("BBRHFSSData::Load", "BBR013", FatalException,
                ("Bad numeric field " + std::string(what) + " = '" + s + "' in " + path).c_str());
    return 0.;
  }
}

// A non-empty row must have exactly n fields: a short row used to be skipped
// silently (BBR013, as the Python mirror reports it).
bool FieldCountOK(const std::vector<std::string>& v, std::size_t n, const G4String& path)
{
  if (v.size() == n) return true;
  G4ExceptionDescription ed;
  ed << "Row with " << v.size() << " fields, expected " << n << ", in " << path;
  G4Exception("BBRHFSSData::Load", "BBR013", FatalException, ed);
  return false;
}

// C1: the CSV header must list exactly the sidecar's columns.
void CheckHeader(const std::string& header, const std::vector<std::string>& want, const G4String& path)
{
  if (SplitCSV(header) != want) {
    G4ExceptionDescription ed;
    ed << "Header of " << path << " differs from the columns its sidecar declares "
          "(BBRHFSSData reads the columns by position).";
    G4Exception("BBRHFSSData::Load", "BBR013", FatalException, ed);
  }
}

bool Within(G4double v, G4double lo, G4double hi, G4double tol) { return v >= lo - tol && v <= hi + tol; }

} // namespace

// ---------------------------------------------------------------------------

BBRHFSSData::BBRHFSSData(const G4String& baseDir, const G4String& dirStem,
                         G4double expectedFreqGHz, const BBRDatasetSidecar* sidecar)
  : fFreqGHz(expectedFreqGHz)
{
  auto dir0 = baseDir + "/" + dirStem + "_Ephi=0";
  auto dir1 = baseDir + "/" + dirStem + "_Ephi=1";
  LoadFarField (dir0 + "/far_field.csv",  0, sidecar);
  LoadFarField (dir1 + "/far_field.csv",  1, sidecar);
  LoadWaveguide(dir0 + "/waveguide.csv",  0, sidecar);
  LoadWaveguide(dir1 + "/waveguide.csv",  1, sidecar);

  if (fData.empty())
    G4Exception("BBRHFSSData", "BBR000", FatalException,
                "No angle datasets loaded — check dataDir path.");

  // Both CSVs use fData[key] (inserting) while loading, so an incidence key
  // present in one file but not the other leaves an empty table that the
  // samplers would index out of bounds. Refuse such a dataset up front.
  for (const auto& [key, ds] : fData) {
    if (ds.farField.empty() || ds.exitPoints.empty()) {
      G4ExceptionDescription ed;
      ed << "Dataset " << dirStem << " key (IWavePhi=" << key.first
         << ", IWaveTheta=" << key.second << ") has " << ds.farField.size()
         << " far-field rows and " << ds.exitPoints.size()
         << " exit-point rows; far_field.csv and waveguide.csv must share the "
            "same incidence keys.";
      G4Exception("BBRHFSSData", "BBR007", FatalException, ed);
    }
  }
  if (sidecar) CheckAgainstSidecar(*sidecar, dirStem);

  // Polarization cross term (Wang eq. 58 applied to eq. 53): the transmitted power of a
  // mixed polarization is |Et E0 + Ep E1|^2 integrated over the exit face, which adds
  // 2 Et Ep Re<E0,E1> to Et^2 T0 + Ep^2 T1. rho is that overlap, normalized.
  //
  // Load-time normalization: the largest transmittance over linear polarizations is
  // the top eigenvalue of [[T0, c], [c, T1]], c = sqrt(T0 T1) Re rho. HFSS port
  // normalization can put it above 1 (raw T0 = 1.0545 at (0,180) of the 52 um gap;
  // T0 + T1 = 1.025 in phase at (45,180)). Such a key has T0 and T1 divided by it, so
  // its maximum is 1 and T1/T0 is kept; a key at or below 1 is left as it is.
  for (auto& [key, ds] : fData) {
    std::complex<G4double> x(0, 0); G4double p0 = 0, p1 = 0;
    for (const auto& e : ds.exitPoints) {
      const std::complex<G4double> a[3] = {{e.Ex_re_0, e.Ex_im_0}, {e.Ey_re_0, e.Ey_im_0}, {e.Ez_re_0, e.Ez_im_0}};
      const std::complex<G4double> b[3] = {{e.Ex_re_1, e.Ex_im_1}, {e.Ey_re_1, e.Ey_im_1}, {e.Ez_re_1, e.Ez_im_1}};
      for (int i = 0; i < 3; ++i) { x += a[i] * std::conj(b[i]); p0 += std::norm(a[i]); p1 += std::norm(b[i]); }
    }
    ds.rho_re = (p0 > 0 && p1 > 0) ? (x / std::sqrt(p0 * p1)).real() : 0.;

    const G4double c = std::sqrt(ds.T_Ephi0 * ds.T_Ephi1) * ds.rho_re;
    const G4double lam = 0.5 * (ds.T_Ephi0 + ds.T_Ephi1) + std::hypot(0.5 * (ds.T_Ephi0 - ds.T_Ephi1), c);
    if (lam > 1.) {
      G4cout << "[BBR] HFSS " << dirStem << " key (" << key.first << ", " << key.second
             << "): max transmittance " << lam
             << " > 1 (port-normalization artefact), normalized to 1" << G4endl;
      ds.T_Ephi0 /= lam;
      ds.T_Ephi1 /= lam;
    }
  }
}

// ---------------------------------------------------------------------------
// Frequency parsing / validation
// ---------------------------------------------------------------------------

// Accepts "<number><MHz|GHz|THz>" (case-insensitive, surrounding blanks ok).
G4double BBRHFSSData::ParseFrequencyGHz(const std::string& tokenIn)
{
  std::string t;
  for (char c : tokenIn)
    if (!std::isspace(static_cast<unsigned char>(c))) t += c;
  if (t.size() < 4) return -1.;

  std::string unit = t.substr(t.size() - 3);
  for (auto& c : unit) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  G4double mult;
  if      (unit == "mhz") mult = 1e-3;
  else if (unit == "ghz") mult = 1.;
  else if (unit == "thz") mult = 1e3;
  else return -1.;

  const std::string num = t.substr(0, t.size() - 3);
  char* end = nullptr;
  const double v = std::strtod(num.c_str(), &end);
  if (num.empty() || end == nullptr || *end != '\0' || !(v > 0.)) return -1.;
  return v * mult;
}

// The frequency a dataset is USED as comes from its directory name; this makes
// sure the file contents agree, so a copied or renamed dataset cannot silently
// stand in for another frequency.
void BBRHFSSData::CheckFrequencyColumn(const G4String& path,
                                       const std::string& token) const
{
  const G4double f = ParseFrequencyGHz(token);
  if (f < 0. || std::abs(f - fFreqGHz) > 1e-3 * fFreqGHz) {
    G4ExceptionDescription ed;
    ed << "Freq column '" << token << "' in " << path
       << " does not match the directory frequency " << fFreqGHz
       << " GHz (tolerance 0.1 %).";
    G4Exception("BBRHFSSData::CheckFrequencyColumn", "BBR009", FatalException, ed);
  }
}

// ---------------------------------------------------------------------------
// CSV loading
// ---------------------------------------------------------------------------

// far_field.csv columns:
//   Freq(0) Ephi(1) IWavePhi(2) IWaveTheta(3) Phi(4) Theta(5)
//   rEphi_real(6) rEphi_imag(7) rEtheta_real(8) rEtheta_imag(9)
void BBRHFSSData::LoadFarField(const G4String& path, int ephi_flag, const BBRDatasetSidecar* sc)
{
  std::ifstream f(path);
  if (!f)
    G4Exception("BBRHFSSData::LoadFarField", "BBR001", FatalException,
                ("Cannot open: " + path).c_str());

  std::string line;
  std::getline(f, line); // skip header
  if (sc) CheckHeader(line, sc->farFieldColumns, path);

  // Row counter per key: the Ephi=1 rows of a key are paired with its Ephi=0
  // FarFieldPoints by position, and each pair must share (Phi, Theta) (BBR012).
  std::map<std::pair<G4double,G4double>, std::size_t> rowIdx;

  while (std::getline(f, line)) {
    auto v = SplitCSV(line);
    if (v.empty()) continue;                    // a blank line (also a bare CRLF)
    if (!FieldCountOK(v, 10, path)) continue;
    CheckFrequencyColumn(path, v[0]);

    G4double iwPhi = Num(v[2], path, "IWavePhi");
    G4double iwThe = Num(v[3], path, "IWaveTheta");
    G4double phi   = Num(v[4], path, "Phi");
    G4double theta = Num(v[5], path, "Theta");
    G4double Epr   = Num(v[6], path, "rEphi_real");
    G4double Epi   = Num(v[7], path, "rEphi_imag");
    G4double Etr   = Num(v[8], path, "rEtheta_real");
    G4double Eti   = Num(v[9], path, "rEtheta_imag");

    auto key = MakeKey(iwPhi, iwThe);
    auto& ds = fData[key];

    if (ephi_flag == 0) {
      FarFieldPoint fp{};
      fp.phi_deg   = phi;  fp.theta_deg   = theta;
      fp.rEphi_re_0 = Epr; fp.rEphi_im_0 = Epi;
      fp.rEtheta_re_0 = Etr; fp.rEtheta_im_0 = Eti;
      ds.farField.push_back(fp);
    } else {
      auto& idx = rowIdx[key];
      if (idx >= ds.farField.size()) {
        G4ExceptionDescription ed;
        ed << "More Ephi=1 than Ephi=0 far-field rows for key (" << iwPhi << ", "
           << iwThe << "): " << path;
        G4Exception("BBRHFSSData::LoadFarField", "BBR012", FatalException, ed);
        return;
      }
      auto& fp = ds.farField[idx++];
      if (std::abs(fp.phi_deg - phi) > 1e-9 || std::abs(fp.theta_deg - theta) > 1e-9) {
        G4ExceptionDescription ed;
        ed << "Ephi=1 far-field row " << (idx - 1) << " of key (" << iwPhi << ", " << iwThe
           << ") is at (Phi, Theta) = (" << phi << ", " << theta
           << ") but the Ephi=0 row is at (" << fp.phi_deg << ", " << fp.theta_deg
           << "): " << path;
        G4Exception("BBRHFSSData::LoadFarField", "BBR012", FatalException, ed);
      }
      fp.rEphi_re_1 = Epr; fp.rEphi_im_1 = Epi;
      fp.rEtheta_re_1 = Etr; fp.rEtheta_im_1 = Eti;
    }
  }

  if (ephi_flag == 1) {
    for (const auto& [key, ds] : fData) {
      const auto it = rowIdx.find(key);
      const std::size_t n1 = (it == rowIdx.end()) ? 0 : it->second;
      if (n1 < ds.farField.size()) {
        G4ExceptionDescription ed;
        ed << "Key (" << key.first << ", " << key.second << ") has " << ds.farField.size()
           << " Ephi=0 but " << n1 << " Ephi=1 far-field rows: " << path;
        G4Exception("BBRHFSSData::LoadFarField", "BBR012", FatalException, ed);
      }
    }
  }
}

// waveguide.csv columns:
//   Freq(0) Ephi(1) IWavePhi(2) IWaveTheta(3) OutgoingPower(4) IngoingPower(5)
//   X(6) Y(7) Z(8) Ex_real(9) Ey_real(10) Ez_real(11) Ex_imag(12) Ey_imag(13) Ez_imag(14)
void BBRHFSSData::LoadWaveguide(const G4String& path, int ephi_flag, const BBRDatasetSidecar* sc)
{
  std::ifstream f(path);
  if (!f)
    G4Exception("BBRHFSSData::LoadWaveguide", "BBR002", FatalException,
                ("Cannot open: " + path).c_str());

  std::string line;
  std::getline(f, line); // skip header
  if (sc) CheckHeader(line, sc->exitFieldColumns, path);

  // Ephi=1 exit points are paired with the Ephi=0 ones by position and must
  // sit at the same (X, Y, Z) (BBR012), as for the far field.
  std::map<std::pair<G4double,G4double>, std::size_t> rowIdx;
  std::map<std::pair<G4double,G4double>, bool>         tSet;

  while (std::getline(f, line)) {
    auto v = SplitCSV(line);
    if (v.empty()) continue;                    // a blank line (also a bare CRLF)
    if (!FieldCountOK(v, 15, path)) continue;
    CheckFrequencyColumn(path, v[0]);

    G4double iwPhi  = Num(v[2], path, "IWavePhi");
    G4double iwThe  = Num(v[3], path, "IWaveTheta");
    G4double outPow = Num(v[4], path, "OutgoingPower");
    G4double inPow  = Num(v[5], path, "IngoingPower");
    G4double x      = Num(v[6], path, "X");
    G4double y      = Num(v[7], path, "Y");
    G4double z      = Num(v[8], path, "Z");
    G4double Ex_re  = Num(v[9], path, "Ex_real");
    G4double Ey_re  = Num(v[10], path, "Ey_real");
    G4double Ez_re  = Num(v[11], path, "Ez_real");
    G4double Ex_im  = Num(v[12], path, "Ex_imag");
    G4double Ey_im  = Num(v[13], path, "Ey_imag");
    G4double Ez_im  = Num(v[14], path, "Ez_imag");

    auto key = MakeKey(iwPhi, iwThe);
    auto& ds = fData[key];

    // Transmittance is constant per (IWavePhi, IWaveTheta): read from first row.
    // The raw ratio can exceed 1 (HFSS port normalization); the constructor
    // normalizes each key once rho is known.
    if (!tSet[key]) {
      const G4double T = (inPow > 0.) ? outPow / inPow : 0.;
      if (ephi_flag == 0) ds.T_Ephi0 = T;
      else                ds.T_Ephi1 = T;
      tSet[key] = true;
    }

    if (ephi_flag == 0) {
      ExitPoint ep{};
      ep.x = x; ep.y = y; ep.z = z;
      ep.Ex_re_0 = Ex_re; ep.Ex_im_0 = Ex_im;
      ep.Ey_re_0 = Ey_re; ep.Ey_im_0 = Ey_im;
      ep.Ez_re_0 = Ez_re; ep.Ez_im_0 = Ez_im;
      ds.exitPoints.push_back(ep);
    } else {
      auto& idx = rowIdx[key];
      if (idx >= ds.exitPoints.size()) {
        G4ExceptionDescription ed;
        ed << "More Ephi=1 than Ephi=0 exit points for key (" << iwPhi << ", "
           << iwThe << "): " << path;
        G4Exception("BBRHFSSData::LoadWaveguide", "BBR012", FatalException, ed);
        return;
      }
      auto& ep = ds.exitPoints[idx++];
      if (std::abs(ep.x - x) > 1e-12 || std::abs(ep.y - y) > 1e-12 || std::abs(ep.z - z) > 1e-12) {
        G4ExceptionDescription ed;
        ed << "Ephi=1 exit point " << (idx - 1) << " of key (" << iwPhi << ", " << iwThe
           << ") is at (" << x << ", " << y << ", " << z << ") but the Ephi=0 point is at ("
           << ep.x << ", " << ep.y << ", " << ep.z << "): " << path;
        G4Exception("BBRHFSSData::LoadWaveguide", "BBR012", FatalException, ed);
      }
      ep.Ex_re_1 = Ex_re; ep.Ex_im_1 = Ex_im;
      ep.Ey_re_1 = Ey_re; ep.Ey_im_1 = Ey_im;
      ep.Ez_re_1 = Ez_re; ep.Ez_im_1 = Ez_im;
    }
  }

  if (ephi_flag == 1) {
    for (const auto& [key, ds] : fData) {
      const auto it = rowIdx.find(key);
      const std::size_t n1 = (it == rowIdx.end()) ? 0 : it->second;
      if (n1 < ds.exitPoints.size()) {
        G4ExceptionDescription ed;
        ed << "Key (" << key.first << ", " << key.second << ") has " << ds.exitPoints.size()
           << " Ephi=0 but " << n1 << " Ephi=1 exit points: " << path;
        G4Exception("BBRHFSSData::LoadWaveguide", "BBR012", FatalException, ed);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Lookup and query
// ---------------------------------------------------------------------------

// C2-C5 (validation/README.md, Dataset sidecars): every key on the declared grid; per key
// the declared far-field and exit row counts; every far-field and exit row in
// the declared ranges; the distinct-value counts; X = 0; and every exit point
// inside the declared cross-section.
void BBRHFSSData::CheckAgainstSidecar(const BBRDatasetSidecar& sc, const G4String& dirStem) const
{
  auto fail = [&](const char* code, const std::string& what) {
    G4ExceptionDescription ed;
    ed << "Dataset " << dirStem << " disagrees with its sidecar " << sc.path << ": " << what;
    G4Exception("BBRHFSSData::CheckAgainstSidecar", code, FatalException, ed);
  };
  auto num = [](G4double v) { std::ostringstream s; s << v; return s.str(); };   // 1e-12, not to_string's 0.000000
  auto fmt = [&](G4double a, G4double b) { return "(" + num(a) + ", " + num(b) + ")"; };
  auto declared = [](G4double v, const std::vector<G4double>& list) {
    for (const G4double d : list) if (RoundDeg(d) == v) return true;
    return false;
  };
  const G4double yTol = 1e-9 * std::max(std::abs(sc.exitY.min), std::abs(sc.exitY.max));
  const G4double zTol = 1e-9 * std::max(std::abs(sc.exitZ.min), std::abs(sc.exitZ.max));
  std::set<G4double> phis, thetas, ys, zs;
  std::size_t outside = 0;
  G4double outY = 0., outZ = 0.;
  for (const auto& [key, ds] : fData) {
    const std::string k = fmt(key.first, key.second);
    if (!declared(key.first, sc.incidentPhiDeg) || !declared(key.second, sc.incidentThetaDeg))
      fail("BBR007", "incidence key " + k + " is not on the declared incident_phi_deg x incident_theta_deg grid");
    if (G4int(ds.farField.size()) != sc.farFieldPointsPerKey)
      fail("BBR012", "key " + k + " has " + std::to_string(ds.farField.size()) + " far-field rows; points_per_key is " +
                     std::to_string(sc.farFieldPointsPerKey));
    if (G4int(ds.exitPoints.size()) != sc.exitPointsPerKey)
      fail("BBR012", "key " + k + " has " + std::to_string(ds.exitPoints.size()) + " exit points; points_per_key_retained is " +
                     std::to_string(sc.exitPointsPerKey));
    for (const auto& fp : ds.farField) {
      if (!Within(fp.phi_deg, sc.farFieldPhi.min, sc.farFieldPhi.max, 1e-9) ||
          !Within(fp.theta_deg, sc.farFieldTheta.min, sc.farFieldTheta.max, 1e-9))
        fail("BBR012", "far-field row (Phi, Theta) = " + fmt(fp.phi_deg, fp.theta_deg) + " of key " + k +
                       " lies outside the declared ranges");
      phis.insert(fp.phi_deg);
      thetas.insert(fp.theta_deg);
    }
    for (const auto& ep : ds.exitPoints) {
      if (std::abs(ep.x) > 1e-12)
        fail("BBR012", "an exit point of key " + k + " has X = " + num(ep.x) + " m; the sidecar declares x_e = 0");
      if (!Within(ep.y, sc.exitY.min, sc.exitY.max, yTol) || !Within(ep.z, sc.exitZ.min, sc.exitZ.max, zTol))
        fail("BBR012", "exit point (Y, Z) = " + fmt(ep.y, ep.z) + " m of key " + k + " lies outside the declared grid");
      ys.insert(ep.y);
      zs.insert(ep.z);
      if (!sc.crossSection.Contains(ep.y, ep.z)) {
        // outside_points "zero": a point with no field in either polarization can never be sampled.
        const bool noField = ep.Ex_re_0 == 0. && ep.Ex_im_0 == 0. && ep.Ey_re_0 == 0. && ep.Ey_im_0 == 0. &&
                             ep.Ez_re_0 == 0. && ep.Ez_im_0 == 0. && ep.Ex_re_1 == 0. && ep.Ex_im_1 == 0. &&
                             ep.Ey_re_1 == 0. && ep.Ey_im_1 == 0. && ep.Ez_re_1 == 0. && ep.Ez_im_1 == 0.;
        if (!(sc.outsidePoints == "zero" && noField)) { ++outside; outY = ep.y; outZ = ep.z; }
      }
    }
  }
  if (G4int(phis.size()) != sc.farFieldPhi.count || G4int(thetas.size()) != sc.farFieldTheta.count)
    fail("BBR012", "the far field has " + std::to_string(phis.size()) + " Phi and " + std::to_string(thetas.size()) +
                   " Theta values; the sidecar declares " + std::to_string(sc.farFieldPhi.count) + " and " +
                   std::to_string(sc.farFieldTheta.count));
  if (G4int(ys.size()) != sc.exitY.count || G4int(zs.size()) != sc.exitZ.count)
    fail("BBR012", "the exit grid has " + std::to_string(ys.size()) + " Y and " + std::to_string(zs.size()) +
                   " Z values; the sidecar declares " + std::to_string(sc.exitY.count) + " and " + std::to_string(sc.exitZ.count));
  if (outside)
    fail("BBR025", std::to_string(outside) + " exit point(s) lie outside the declared cross-section, e.g. (Y, Z) = " +
                   fmt(outY, outZ) + " m");
}

const BBRHFSSData::AngleDataset& BBRHFSSData::FindDataset(
    G4double iwavePhi_deg, G4double iwaveTheta_deg) const
{
  const AngleDataset* best = nullptr;
  G4double bestDist2 = 1e30;
  for (const auto& [key, ds] : fData) {
    G4double dp = key.first  - iwavePhi_deg;
    G4double dt = key.second - iwaveTheta_deg;
    G4double d2 = dp*dp + dt*dt;
    if (d2 < bestDist2) { bestDist2 = d2; best = &ds; }
  }
  // fData is non-empty by construction (BBR000); best is null only for NaN angles.
  return best ? *best : fData.begin()->second;
}

G4double BBRHFSSData::GetTransmittance(G4double E_theta, G4double E_phi,
                                       G4double iwavePhi_deg,
                                       G4double iwaveTheta_deg) const
{
  const auto& ds = FindDataset(iwavePhi_deg, iwaveTheta_deg);
  // The load-time normalization keeps the largest T over linear polarizations at
  // or below 1; the clamp guards against rounding and non-unit (E_theta, E_phi).
  G4double T = E_theta*E_theta * ds.T_Ephi0 + E_phi*E_phi * ds.T_Ephi1
             + 2.*E_theta*E_phi * std::sqrt(ds.T_Ephi0 * ds.T_Ephi1) * ds.rho_re;
  return std::min(1., std::max(0., T));
}

// ---------------------------------------------------------------------------

G4ThreeVector BBRHFSSData::SampleOutgoingDirection(
    G4double E_theta, G4double E_phi,
    G4double iwavePhi_deg, G4double iwaveTheta_deg,
    const G4ThreeVector& phi_hat,
    const G4ThreeVector& theta_hat,
    const G4ThreeVector& normal_hat,
    G4ThreeVector& pol_out) const
{
  const auto& ds = FindDataset(iwavePhi_deg, iwaveTheta_deg);
  const auto& ff = ds.farField;
  const std::size_t N = ff.size();

  // Runtime CDF: combined power |E_theta·F₀ + E_phi·F₁|² per far-field point,
  // weighted by sinT to account for solid angle dΩ = sinT·dT·dPhi (this also
  // zeroes the Theta = 0 / 180 rows, which all map to the same direction).
  // Directions lying in the exit-face plane (normal component below
  // kMinNormalComponent) get zero weight: see the constant's comment.
  std::vector<G4double> cdf(N);
  G4double sum = 0.;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& fp = ff[i];
    G4double sinT   = std::sin(fp.theta_deg * CLHEP::pi / 180.);
    G4double cosN   = sinT * std::cos(fp.phi_deg * CLHEP::pi / 180.);
    if (cosN < kMinNormalComponent) { cdf[i] = sum; continue; }
    G4double Eth_re = E_theta*fp.rEtheta_re_0 + E_phi*fp.rEtheta_re_1;
    G4double Eth_im = E_theta*fp.rEtheta_im_0 + E_phi*fp.rEtheta_im_1;
    G4double Eph_re = E_theta*fp.rEphi_re_0   + E_phi*fp.rEphi_re_1;
    G4double Eph_im = E_theta*fp.rEphi_im_0   + E_phi*fp.rEphi_im_1;
    sum    += sinT * (Eth_re*Eth_re + Eth_im*Eth_im + Eph_re*Eph_re + Eph_im*Eph_im);
    cdf[i]  = sum;
  }

  std::size_t k;
  if (sum > 0.) {
    for (auto& c : cdf) c /= sum;
    G4double U = G4UniformRand();
    k = (std::size_t)(
        std::lower_bound(cdf.begin(), cdf.end(), U) - cdf.begin());
    if (k >= N) k = N - 1;
  } else {
    // Degenerate dataset (all amplitudes zero): fall back to a uniform pick.
    k = std::min<std::size_t>(N - 1, (std::size_t)(G4UniformRand() * N));
  }

  const auto& fp = ff[k];
  G4double T_rad = fp.theta_deg * CLHEP::pi / 180.;
  G4double P_rad = fp.phi_deg   * CLHEP::pi / 180.;
  G4double sinT = std::sin(T_rad), cosT = std::cos(T_rad);
  G4double sinP = std::sin(P_rad), cosP = std::cos(P_rad);

  G4ThreeVector dir_out = sinT*cosP*normal_hat + sinT*sinP*theta_hat + cosT*phi_hat;

  // Outgoing spherical basis vectors at (T,P).
  G4ThreeVector eTh_out =  cosT*cosP*normal_hat + cosT*sinP*theta_hat - sinT*phi_hat;
  G4ThreeVector ePh_out = -sinP*normal_hat       + cosP*theta_hat;

  // Geant4 wants a real polarization vector, but the sampled far-field
  // amplitude is complex (elliptical in general). Use the major axis of the
  // polarization ellipse: |Re[(a êθ + b êφ) e^{iψ}]| is maximal at
  // ψ = −arg(a² + b²)/2. (Using only the real parts — the previous behaviour
  // — fails when the amplitudes are predominantly imaginary.)
  std::complex<G4double> a(E_theta*fp.rEtheta_re_0 + E_phi*fp.rEtheta_re_1,
                           E_theta*fp.rEtheta_im_0 + E_phi*fp.rEtheta_im_1);
  std::complex<G4double> b(E_theta*fp.rEphi_re_0   + E_phi*fp.rEphi_re_1,
                           E_theta*fp.rEphi_im_0   + E_phi*fp.rEphi_im_1);
  std::complex<G4double> s = a*a + b*b;
  G4double psi = (std::abs(s) > 0.) ? -0.5 * std::arg(s) : 0.;
  std::complex<G4double> rot = std::polar(1., psi);
  pol_out = std::real(a*rot)*eTh_out + std::real(b*rot)*ePh_out;
  if (pol_out.mag() > 1e-30) pol_out = pol_out.unit();
  else                        pol_out = eTh_out;  // degenerate fallback

  return dir_out.unit();
}

// ---------------------------------------------------------------------------

G4ThreeVector BBRHFSSData::SampleExitPosition(
    G4double E_theta, G4double E_phi,
    G4double iwavePhi_deg, G4double iwaveTheta_deg,
    const G4ThreeVector& exit_face_center,
    const G4ThreeVector& crack_x,
    const G4ThreeVector& crack_y) const
{
  const auto& ds  = FindDataset(iwavePhi_deg, iwaveTheta_deg);
  const auto& eps = ds.exitPoints;
  const std::size_t M = eps.size();

  // Runtime CDF: combined |E_theta·E₀(x,y) + E_phi·E₁(x,y)|² per exit point.
  std::vector<G4double> cdf(M);
  G4double sum = 0.;
  for (std::size_t j = 0; j < M; ++j) {
    const auto& ep = eps[j];
    G4double px_re = E_theta*ep.Ex_re_0 + E_phi*ep.Ex_re_1;
    G4double px_im = E_theta*ep.Ex_im_0 + E_phi*ep.Ex_im_1;
    G4double py_re = E_theta*ep.Ey_re_0 + E_phi*ep.Ey_re_1;
    G4double py_im = E_theta*ep.Ey_im_0 + E_phi*ep.Ey_im_1;
    G4double pz_re = E_theta*ep.Ez_re_0 + E_phi*ep.Ez_re_1;
    G4double pz_im = E_theta*ep.Ez_im_0 + E_phi*ep.Ez_im_1;
    sum    += px_re*px_re + px_im*px_im + py_re*py_re + py_im*py_im
            + pz_re*pz_re + pz_im*pz_im;
    cdf[j]  = sum;
  }

  std::size_t j;
  if (sum > 0.) {
    for (auto& c : cdf) c /= sum;
    G4double U = G4UniformRand();
    j = (std::size_t)(
        std::lower_bound(cdf.begin(), cdf.end(), U) - cdf.begin());
    if (j >= M) j = M - 1;
  } else {
    j = std::min<std::size_t>(M - 1, (std::size_t)(G4UniformRand() * M));
  }

  const auto& ep = eps[j];
  // HFSS coordinates are SI (meters); Geant4's base unit is mm, so multiply by
  // CLHEP::m. The CSV points are in the exit coordinate system (x_e, y_e, z_e) =
  // (Z, Y, -X)_HFSS = (normal, theta_hat, phi_hat): the exit plane is x_e = 0 (CSV X
  // column = 0), CSV Y runs along crack_x (theta_f) and CSV Z along crack_y (phi_f).
  return exit_face_center + (ep.y * CLHEP::m) * crack_x
                          + (ep.z * CLHEP::m) * crack_y;
}
