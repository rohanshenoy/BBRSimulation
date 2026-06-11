#include "BBRHFSSData.hh"

#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> SplitCSV(const std::string& line)
{
  std::vector<std::string> v;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, ',')) v.push_back(tok);
  return v;
}

// Stable floating-point map key: round to nearest 0.01 degree.
G4double RoundDeg(G4double d) { return std::round(d * 100.) / 100.; }

std::pair<G4double, G4double> MakeKey(G4double phi, G4double theta)
{
  return {RoundDeg(phi), RoundDeg(theta)};
}

} // namespace

// ---------------------------------------------------------------------------

BBRHFSSData::BBRHFSSData(const G4String& baseDir, const G4String& datasetId)
{
  auto dir0 = baseDir + "/" + datasetId + "_Ephi=0";
  auto dir1 = baseDir + "/" + datasetId + "_Ephi=1";
  LoadFarField (dir0 + "/far_field.csv",  0);
  LoadFarField (dir1 + "/far_field.csv",  1);
  LoadWaveguide(dir0 + "/waveguide.csv",  0);
  LoadWaveguide(dir1 + "/waveguide.csv",  1);

  if (fData.empty())
    G4Exception("BBRHFSSData", "BBR000", FatalException,
                "No angle datasets loaded — check dataDir path.");
}

// ---------------------------------------------------------------------------
// CSV loading
// ---------------------------------------------------------------------------

// far_field.csv columns:
//   Freq(0) Ephi(1) IWavePhi(2) IWaveTheta(3) Phi(4) Theta(5)
//   rEphi_real(6) rEphi_imag(7) rEtheta_real(8) rEtheta_imag(9)
void BBRHFSSData::LoadFarField(const G4String& path, int ephi_flag)
{
  std::ifstream f(path);
  if (!f)
    G4Exception("BBRHFSSData::LoadFarField", "BBR001", FatalException,
                ("Cannot open: " + path).c_str());

  std::string line;
  std::getline(f, line); // skip header

  // Row counter per key, used to match ephi=1 rows to the ephi=0 FarFieldPoints
  // by insertion order (both CSVs have the same angular sweep order).
  std::map<std::pair<G4double,G4double>, std::size_t> rowIdx;

  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto v = SplitCSV(line);
    if (v.size() < 10) continue;

    G4double iwPhi = std::stod(v[2]);
    G4double iwThe = std::stod(v[3]);
    G4double phi   = std::stod(v[4]);
    G4double theta = std::stod(v[5]);
    G4double Epr   = std::stod(v[6]);
    G4double Epi   = std::stod(v[7]);
    G4double Etr   = std::stod(v[8]);
    G4double Eti   = std::stod(v[9]);

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
      if (idx < ds.farField.size()) {
        auto& fp = ds.farField[idx++];
        fp.rEphi_re_1 = Epr; fp.rEphi_im_1 = Epi;
        fp.rEtheta_re_1 = Etr; fp.rEtheta_im_1 = Eti;
      }
    }
  }
}

// waveguide.csv columns:
//   Freq(0) Ephi(1) IWavePhi(2) IWaveTheta(3) OutgoingPower(4) IngoingPower(5)
//   X(6) Y(7) Z(8) Ex_real(9) Ey_real(10) Ez_real(11) Ex_imag(12) Ey_imag(13) Ez_imag(14)
void BBRHFSSData::LoadWaveguide(const G4String& path, int ephi_flag)
{
  std::ifstream f(path);
  if (!f)
    G4Exception("BBRHFSSData::LoadWaveguide", "BBR002", FatalException,
                ("Cannot open: " + path).c_str());

  std::string line;
  std::getline(f, line); // skip header

  std::map<std::pair<G4double,G4double>, std::size_t> rowIdx;
  std::map<std::pair<G4double,G4double>, bool>         tSet;

  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto v = SplitCSV(line);
    if (v.size() < 15) continue;

    G4double iwPhi  = std::stod(v[2]);
    G4double iwThe  = std::stod(v[3]);
    G4double outPow = std::stod(v[4]);
    G4double inPow  = std::stod(v[5]);
    G4double x      = std::stod(v[6]);
    G4double y      = std::stod(v[7]);
    G4double z      = std::stod(v[8]);
    G4double Ex_re  = std::stod(v[9]);
    G4double Ey_re  = std::stod(v[10]);
    G4double Ez_re  = std::stod(v[11]);
    G4double Ex_im  = std::stod(v[12]);
    G4double Ey_im  = std::stod(v[13]);
    G4double Ez_im  = std::stod(v[14]);

    auto key = MakeKey(iwPhi, iwThe);
    auto& ds = fData[key];

    // Transmittance is constant per (IWavePhi, IWaveTheta): read from first row.
    if (!tSet[key]) {
      G4double T = (inPow > 0.) ? outPow / inPow : 0.;
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
      if (idx < ds.exitPoints.size()) {
        auto& ep = ds.exitPoints[idx++];
        ep.Ex_re_1 = Ex_re; ep.Ex_im_1 = Ex_im;
        ep.Ey_re_1 = Ey_re; ep.Ey_im_1 = Ey_im;
        ep.Ez_re_1 = Ez_re; ep.Ez_im_1 = Ez_im;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Lookup and query
// ---------------------------------------------------------------------------

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
  if (!best)
    G4Exception("BBRHFSSData::FindDataset", "BBR003", FatalException, "fData is empty");
  return *best;
}

G4double BBRHFSSData::GetTransmittance(G4double E_theta, G4double E_phi,
                                       G4double iwavePhi_deg,
                                       G4double iwaveTheta_deg) const
{
  const auto& ds = FindDataset(iwavePhi_deg, iwaveTheta_deg);
  // HFSS power ratios can slightly exceed 1 (e.g. 1.055 at normal incidence,
  // a numerical artefact of the port normalization) — clamp to a probability.
  G4double T = E_theta*E_theta * ds.T_Ephi0 + E_phi*E_phi * ds.T_Ephi1;
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
  // weighted by sinT to account for solid angle dΩ = sinT·dT·dPhi.
  // that all map to the same physical direction.
  std::vector<G4double> cdf(N);
  G4double sum = 0.;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& fp = ff[i];
    G4double sinT   = std::sin(fp.theta_deg * CLHEP::pi / 180.);
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
  // HFSS coordinates are SI (meters). Geant4 base unit is mm → multiply by CLHEP::m.
  // HFSS model: x=b(gap), y=long, z=propagation. Exit face (X=0 in CSV):
  //   CSV Y → y_model → theta_hat (world ŷ)   via crack_x
  //   CSV Z → x_model → phi_hat   (world ẑ)   via crack_y
  return exit_face_center + (ep.y * CLHEP::m) * crack_x
                          + (ep.z * CLHEP::m) * crack_y;
}
