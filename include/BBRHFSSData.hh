#ifndef BBRHFSSData_hh
#define BBRHFSSData_hh

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include <map>
#include <utility>
#include <vector>

// Loads HFSS far-field and waveguide CSVs and samples the diffraction boundary
// condition per Wang (2025) eqs. 54-58.
//
// Callers pass:
//   E_theta = phat · ê_θ_in  (real scalar, incoming HFSS θ-basis component)
//   E_phi   = phat · ê_φ_in  (real scalar, incoming HFSS φ-basis component)
// where ê_θ_in, ê_φ_in are computed in world frame from (IWaveTheta, IWavePhi)
// by HandleDiffractionBoundary. E_theta²+E_phi²=1 since phat⊥khat.
//
// CDFs are built at runtime (not precomputed) so the E_theta×E_phi cross term
// in |E_theta·F₀ + E_phi·F₁|² is handled exactly.
class BBRHFSSData
{
 public:
  // baseDir: path to HFSSSimData directory (e.g. "../HFSSSimData")
  // datasetId: subdirectory prefix, e.g. "InfParallelPlate_crack1Rohan_500GHz"
  //   Loads baseDir/datasetId_Ephi=0/ and baseDir/datasetId_Ephi=1/
  BBRHFSSData(const G4String& baseDir, const G4String& datasetId);

  // Wang eq. 54: T = E_theta²·T₀ + E_phi²·T₁
  G4double GetTransmittance(G4double E_theta, G4double E_phi,
                            G4double iwavePhi_deg, G4double iwaveTheta_deg) const;

  // Wang eq. 56-57: sample outgoing direction; pol_out set as output param.
  G4ThreeVector SampleOutgoingDirection(G4double E_theta, G4double E_phi,
                                        G4double iwavePhi_deg, G4double iwaveTheta_deg,
                                        const G4ThreeVector& phi_hat,
                                        const G4ThreeVector& theta_hat,
                                        const G4ThreeVector& normal_hat,
                                        G4ThreeVector& pol_out) const;

  // Wang eq. 55: sample exit-face position (world frame).
  //   exit_face_center: world position of exit-face center
  //   crack_x: world unit vector for HFSS waveguide Y axis (= theta_hat for standard geometry)
  //   crack_y: world unit vector for HFSS waveguide Z axis (= phi_hat  for standard geometry)
  G4ThreeVector SampleExitPosition(G4double E_theta, G4double E_phi,
                                   G4double iwavePhi_deg, G4double iwaveTheta_deg,
                                   const G4ThreeVector& exit_face_center,
                                   const G4ThreeVector& crack_x,
                                   const G4ThreeVector& crack_y) const;

 private:
  struct FarFieldPoint {
    G4double phi_deg, theta_deg;
    // Complex amplitudes for both HFSS basis datasets (suffix _0 = Ephi=0, _1 = Ephi=1).
    G4double rEtheta_re_0, rEtheta_im_0, rEphi_re_0, rEphi_im_0;
    G4double rEtheta_re_1, rEtheta_im_1, rEphi_re_1, rEphi_im_1;
  };

  struct ExitPoint {
    G4double x, y, z;           // HFSS SI coordinates (meters)
    G4double Ex_re_0, Ex_im_0, Ey_re_0, Ey_im_0, Ez_re_0, Ez_im_0;  // Ephi=0
    G4double Ex_re_1, Ex_im_1, Ey_re_1, Ey_im_1, Ez_re_1, Ez_im_1;  // Ephi=1
  };

  struct AngleDataset {
    std::vector<FarFieldPoint> farField;
    G4double T_Ephi0 = 0.;   // transmittance, θ-polarised (Ephi=0) input
    G4double T_Ephi1 = 0.;   // transmittance, φ-polarised (Ephi=1) input
    std::vector<ExitPoint> exitPoints;
  };

  // Keyed by (RoundedIWavePhi_deg, RoundedIWaveTheta_deg); populated dynamically from CSV.
  std::map<std::pair<G4double, G4double>, AngleDataset> fData;

  // Nearest-neighbour lookup by L2 distance in (phi, theta) degree space.
  const AngleDataset& FindDataset(G4double iwavePhi_deg, G4double iwaveTheta_deg) const;

  // ephi_flag: 0 = Ephi=0 CSV (fill _0 fields), 1 = Ephi=1 CSV (fill _1 fields).
  void LoadFarField(const G4String& path, int ephi_flag);
  void LoadWaveguide(const G4String& path, int ephi_flag);
};

#endif
