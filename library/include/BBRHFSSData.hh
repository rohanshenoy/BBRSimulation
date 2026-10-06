#ifndef BBRHFSSData_hh
#define BBRHFSSData_hh

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#include <map>
#include <utility>
#include <vector>

struct BBRDatasetSidecar;

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
//
// Fatal G4Exception codes raised while loading:
//   BBR000  no incidence key in either CSV
//   BBR001  far_field.csv cannot be opened
//   BBR002  waveguide.csv cannot be opened
//   BBR007  far_field.csv and waveguide.csv have different incidence keys; with
//           a sidecar, also C2 (validation/README.md, Dataset sidecars)
//   BBR009  a row's Freq disagrees with the directory frequency (every row)
//   BBR012  an Ephi=1 row does not match the Ephi=0 row at the same index
//           (X/Y/Z of an exit point, Phi/Theta of a far-field row), or the
//           two files have different row counts for a key; with a sidecar,
//           also C3-C4 (validation/README.md, Dataset sidecars)
//   BBR013  a numeric field is not a number or not finite, or a non-empty
//           row has the wrong number of fields; with a sidecar, also C1
//           (validation/README.md, Dataset sidecars)
//   BBR025  (with a sidecar) an exit point outside the declared cross-section
class BBRHFSSData
{
 public:
  // baseDir: the waveguide data directory (<data root>/waveguides).
  // dirStem: directory name without the "_Ephi=N" suffix, i.e.
  //   "<id>_<freq>GHz" (e.g. "InfParallelPlate_crack1Rohan_500GHz").
  //   Loads baseDir/dirStem_Ephi=0/ and baseDir/dirStem_Ephi=1/.
  // expectedFreqGHz: the frequency parsed from the directory name. Each CSV's
  //   Freq column must agree with it to 0.1 % (BBR009), so a mislabelled or
  //   mis-copied dataset cannot be used under the wrong frequency.
  // sidecar (optional): the dataset's parsed <dirStem>.dataset.json. When
  // given, the loader also runs the sidecar checks C1-C5: header = declared
  // columns (BBR013), keys on the declared grid (BBR007), far-field and
  // exit-grid rows per key, ranges and distinct-value counts (equal to the
  // declared ones; with outside_points "omitted" the exit counts may be
  // smaller, never larger, because the producer declares the whole export
  // lattice and drops the points outside the section), X = 0 (BBR012), every
  // exit point inside the declared cross-section (BBR025). BBRCrackLibrary always
  // passes it; a direct construction without one keeps the pre-sidecar checks.
  BBRHFSSData(const G4String& baseDir, const G4String& dirStem,
              G4double expectedFreqGHz, const BBRDatasetSidecar* sidecar = nullptr);

  // Frequency this dataset was loaded as [GHz].
  G4double GetFrequencyGHz() const { return fFreqGHz; }

  // "500GHz" -> 500, "0.5THz" -> 500, "500000MHz" -> 500; -1 if unparseable.
  // Shared with BBRCrackLibrary, which parses the same token from directory names.
  static G4double ParseFrequencyGHz(const std::string& token);

  // T = E_theta²·T₀ + E_phi²·T₁ + 2·E_theta·E_phi·√(T₀T₁)·Re ρ, clamped to [0, 1]
  // (Wang eq. 58 substituted into eq. 53), with ρ the normalized overlap
  // Σ E₀·E₁* / √(Σ|E₀|² Σ|E₁|²) of the two basis exit fields: the transmitted
  // power is quadratic in the total exit field, so the two powers add without a
  // cross term only when the basis fields are orthogonal over the exit face.
  // T₀ and T₁ are OutgoingPower/IngoingPower of the key's first row; where the
  // largest T over linear polarizations (the top eigenvalue of [[T₀, c], [c, T₁]],
  // c = √(T₀T₁)·Re ρ) exceeds 1, an HFSS port-normalization artefact, the load
  // divides both by it and logs "[BBR] HFSS … normalized to 1".
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
    G4double rho_re = 0.;  // Re of the normalized exit-field overlap <E0,E1>, for the polarization cross term
    std::vector<ExitPoint> exitPoints;
  };

  // Keyed by (RoundedIWavePhi_deg, RoundedIWaveTheta_deg); populated dynamically from CSV.
  std::map<std::pair<G4double, G4double>, AngleDataset> fData;

  G4double fFreqGHz = -1.;   // frequency of this dataset [GHz]

  // Fatal BBR009 if the CSV's Freq column disagrees with fFreqGHz.
  void CheckFrequencyColumn(const G4String& path, const std::string& token) const;

  // Nearest-neighbour lookup by L2 distance in (phi, theta) degree space.
  const AngleDataset& FindDataset(G4double iwavePhi_deg, G4double iwaveTheta_deg) const;

  void CheckAgainstSidecar(const BBRDatasetSidecar& sc, const G4String& dirStem) const;

  // ephi_flag: 0 = Ephi=0 CSV (fill _0 fields), 1 = Ephi=1 CSV (fill _1 fields).
  void LoadFarField(const G4String& path, int ephi_flag, const BBRDatasetSidecar* sc);
  void LoadWaveguide(const G4String& path, int ephi_flag, const BBRDatasetSidecar* sc);
};

#endif
