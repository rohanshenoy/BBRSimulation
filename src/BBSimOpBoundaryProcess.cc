#include "BBSimOpBoundaryProcess.hh"

#include "BBRCrackLibrary.hh"
#include "G4AffineTransform.hh"
#include "G4NavigationHistory.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"

#include <cmath>

BBSimOpBoundaryProcess::BBSimOpBoundaryProcess(const G4String& name)
  : G4WrapperProcess(name)
{}

// PostStepDoIt — intercept TEM_waveguide volumes; fall through otherwise.
// The pass-through path preserves bit-identical output for all other geometry
// (verify.mac regression must still pass).
G4VParticleChange* BBSimOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                                        const G4Step& aStep)
{
  const G4Material* mat2 = aStep.GetPostStepPoint()->GetMaterial();
  if (mat2 && mat2->GetName() == "vacuum_wg")
    return HandleDiffractionBoundary(aTrack, aStep);
  if (mat2) {
    G4MaterialPropertiesTable* mpt = mat2->GetMaterialPropertiesTable();
    if (mpt && mpt->GetProperty("REFLECTIVITY"))
      return HandleReflectanceBoundary(aTrack, aStep);
  }
  return pRegProcess->PostStepDoIt(aTrack, aStep);
}

// ---------------------------------------------------------------------------

G4VParticleChange* BBSimOpBoundaryProcess::HandleDiffractionBoundary(
    const G4Track& aTrack, const G4Step& aStep)
{
  fParticleChange.Initialize(aTrack);

  const G4ThreeVector khat = aTrack.GetMomentumDirection();
  const G4ThreeVector phat = aTrack.GetPolarization();

  // Extract crack-local axes from the volume's rotation in world frame.
  // Convention: local +X = propagation (normal_hat), +Y = long dim (theta_hat),
  //             +Z = gap (phi_hat). InverseTransformAxis maps local → world.
  const G4VTouchable*      touch = aStep.GetPostStepPoint()->GetTouchable();
  const G4AffineTransform& xf    = touch->GetHistory()->GetTopTransform();
  G4ThreeVector normal_hat = xf.InverseTransformAxis(G4ThreeVector(1., 0., 0.));
  G4ThreeVector theta_hat  = xf.InverseTransformAxis(G4ThreeVector(0., 1., 0.));
  G4ThreeVector phi_hat    = xf.InverseTransformAxis(G4ThreeVector(0., 0., 1.));
  // Flip normal_hat to point outward (same half-space as khat, i.e. toward exit face).
  if (khat.dot(normal_hat) < 0.) normal_hat = -normal_hat;

  // Lazy-load the HFSS dataset for this crack volume.
  const G4String volName   = touch->GetVolume()->GetName();
  const G4String datasetId = volName.substr(0, volName.find(':'));
  BBRHFSSData& hfss = BBRCrackLibrary::Instance().Lookup(datasetId);

  // --- incoming angles in crack-local frame (folded into HFSS quarter-symmetry) ---
  // HFSS convention: ẑ_i points OUT of the crack (= normal_hat = +x_world).
  // A photon propagating along +x_world (khat = normal_hat) approaches the
  // exit face from inside, i.e. its k-vector is anti-parallel to ẑ_i →
  // IWaveTheta = 180° (confirmed: T≈1.055 at IWaveTheta=180°, T≈0 at 0°).
  G4double cosVal = std::min(1., std::max(-1., -khat.dot(normal_hat)));
  G4double iwaveTheta_deg = std::acos(cosVal) * (180. / CLHEP::pi);
  // IWavePhi: azimuth in HFSS x̂_i=phi_hat, ŷ_i=-theta_hat plane.
  G4double iwavePhi_raw   = std::atan2(-khat.dot(theta_hat), khat.dot(phi_hat))
                            * (180. / CLHEP::pi);
  G4double iwavePhi_deg   = std::abs(iwavePhi_raw);
  if (iwavePhi_deg > 90.) iwavePhi_deg = 180. - iwavePhi_deg;

  // --- decompose incoming polarization onto HFSS incoming spherical basis ---
  // Standard spherical basis at (IWaveTheta, IWavePhi) in HFSS frame
  // (ẑ_i=normal_hat, x̂_i=phi_hat, ŷ_i=-theta_hat):
  //   ê_θ = -sin(T)*ẑ_i + cos(T)*cos(P)*x̂_i - cos(T)*sin(P)*ŷ_i
  //   ê_φ =            -sin(P)*x̂_i           - cos(P)*ŷ_i
  G4double th = iwaveTheta_deg * (CLHEP::pi / 180.);
  G4double ph = iwavePhi_raw   * (CLHEP::pi / 180.);
  G4ThreeVector eTheta_in = +std::sin(th) * normal_hat
                            - std::cos(th) * std::sin(ph) * theta_hat
                            + std::cos(th) * std::cos(ph) * phi_hat;
  G4ThreeVector ePhi_in   = -std::cos(ph) * theta_hat
                            - std::sin(ph) * phi_hat;

  G4double E_theta = phat.dot(eTheta_in);
  G4double E_phi   = phat.dot(ePhi_in);
  G4double norm    = std::sqrt(E_theta*E_theta + E_phi*E_phi);
  if (norm > 1e-9) { E_theta /= norm; E_phi /= norm; }
  else             { E_theta = M_SQRT1_2; E_phi = M_SQRT1_2; }

  // --- transmittance decision (Wang eq. 54) ---
  G4double T = hfss.GetTransmittance(E_theta, E_phi,
                                     iwavePhi_deg, iwaveTheta_deg);

  const G4bool transmitted = (G4UniformRand() < T);
  ++fNDiffraction;
  if (transmitted) ++fNDiffractionTransmit;
  if (fNDiffraction % 100 == 0)
    G4cout << "[BBR] diffraction events=" << fNDiffraction
           << " T_obs=" << G4double(fNDiffractionTransmit)/fNDiffraction << G4endl;

  G4ThreeVector dir_final;
  G4ThreeVector pos_csv;
  G4bool hasPos = false;

  if (transmitted) {
    // Transmit: sample outgoing direction + polarization (Wang eqs. 56-57).
    G4ThreeVector pol_out;
    G4ThreeVector dir_out = hfss.SampleOutgoingDirection(
        E_theta, E_phi, iwavePhi_deg, iwaveTheta_deg,
        phi_hat, theta_hat, normal_hat, pol_out);

    // Sample exit-face position (Wang eq. 55).
    // Exit face center: traverse the crack slab from the entry point.
    const G4ThreeVector& entryGl = aStep.GetPostStepPoint()->GetPosition();
    G4ThreeVector localEntry     = xf.TransformPoint(entryGl);
    G4ThreeVector localNormal    = xf.TransformAxis(normal_hat);
    // Nudge slightly inside so DistanceToOut is well-defined at the boundary.
    G4ThreeVector localStart     = localEntry + 1e-7*CLHEP::mm * localNormal;
    G4double thickness           = touch->GetSolid()->DistanceToOut(localStart, localNormal);
    G4ThreeVector exitCenter     = xf.InverseTransformPoint(localStart + thickness * localNormal);

    // HFSS waveguide CSV: Y=long dim → +theta_hat, Z=gap/b → phi_hat.
    // (ŷ_i = −ŷ_world sign flip is for angle conventions only, not exit positions.)
    G4ThreeVector pos_out = hfss.SampleExitPosition(
        E_theta, E_phi, iwavePhi_deg, iwaveTheta_deg,
        exitCenter, theta_hat, phi_hat);

    dir_final = dir_out;
    pos_csv   = pos_out;
    hasPos    = true;

    fParticleChange.ProposeMomentumDirection(dir_out);
    fParticleChange.ProposePolarization(pol_out);
    fParticleChange.ProposePosition(pos_out);
    fParticleChange.ProposeTrackStatus(fAlive);
  } else {
    // Reflect: specular flip of momentum and polarization about the crack normal.
    G4ThreeVector dir_ref = (khat - 2.*khat.dot(normal_hat)*normal_hat).unit();
    G4ThreeVector pol_ref = phat  - 2.*phat.dot(normal_hat)*normal_hat;
    if (pol_ref.mag() > 1e-30) pol_ref = pol_ref.unit();
    else                        pol_ref = phi_hat;

    dir_final = dir_ref;

    fParticleChange.ProposeMomentumDirection(dir_ref);
    fParticleChange.ProposePolarization(pol_ref);
    fParticleChange.ProposeTrackStatus(fAlive);
  }

  return &fParticleChange;
}

// ---------------------------------------------------------------------------

G4VParticleChange* BBSimOpBoundaryProcess::HandleReflectanceBoundary(
    const G4Track& aTrack, const G4Step& aStep)
{
  fParticleChange.Initialize(aTrack);

  // Read REFLECTIVITY from Material2's MPT (YYC pattern, adapted for G4 11.4).
  const G4Material*          mat2 = aStep.GetPostStepPoint()->GetMaterial();
  G4MaterialPropertiesTable* mpt  = mat2->GetMaterialPropertiesTable();
  G4MaterialPropertyVector*  rvec =
      mpt->GetProperty("REFLECTIVITY");
  G4double E = aTrack.GetKineticEnergy();   // YYC: thePhotonMomentum
  G4double R = rvec->Value(E);              // YYC: PropertyPointer->Value(thePhotonMomentum)

  // Surface normal pointing away from mat2 back into the vacuum (opposes k).
  // YYC: theGlobalNormal (computed internally by G4OpBoundaryProcess).
  const G4VTouchable*      postTouch = aStep.GetPostStepPoint()->GetTouchable();
  const G4AffineTransform& postXF   = postTouch->GetHistory()->GetTopTransform();
  G4ThreeVector posLocal  = postXF.TransformPoint(
                                aStep.GetPostStepPoint()->GetPosition());
  G4ThreeVector normLocal = postTouch->GetSolid()->SurfaceNormal(posLocal);
  G4ThreeVector nhat      = postXF.InverseTransformAxis(normLocal);
  if (nhat.dot(aTrack.GetMomentumDirection()) > 0.) nhat = -nhat;

  ++fNReflectance;

  // YYC dispatch: rand > R → DoAbsorption(); rand <= R → DoReflection() specular.
  G4double rand = G4UniformRand();
  if (rand > R) {
    ++fNReflectanceAbsorb;
    fParticleChange.ProposeLocalEnergyDeposit(E);
    fParticleChange.ProposeTrackStatus(fStopAndKill);
  } else {
    // k_ref = k − 2(k·n)n;  p_ref = p − 2(p·n)n
    const G4ThreeVector& k = aTrack.GetMomentumDirection(); // YYC: OldMomentum
    const G4ThreeVector& p = aTrack.GetPolarization();      // YYC: OldPolarization
    G4ThreeVector k_ref = (k - 2.*k.dot(nhat)*nhat).unit();
    G4ThreeVector p_ref =  p - 2.*p.dot(nhat)*nhat;
    if (p_ref.mag() > 1e-30) p_ref = p_ref.unit();
    else                      p_ref = k_ref.cross(nhat).unit();
    fParticleChange.ProposeMomentumDirection(k_ref);
    fParticleChange.ProposePolarization(p_ref);
    fParticleChange.ProposeTrackStatus(fAlive);
  }

  if (fNReflectance % 1000 == 0)
    G4cout << "[BBR] reflectance mat=" << mat2->GetName()
           << " N=" << fNReflectance
           << " A_obs=" << G4double(fNReflectanceAbsorb)/fNReflectance
           << " R_theory=" << R << G4endl;

  return &fParticleChange;
}

// ---------------------------------------------------------------------------

G4OpBoundaryProcessStatus BBSimOpBoundaryProcess::GetStatus() const
{
  return GetWrappedProcess()->GetStatus();
}

void BBSimOpBoundaryProcess::SetInvokeSD(G4bool flag)
{
  GetWrappedProcess()->SetInvokeSD(flag);
}

G4OpBoundaryProcess* BBSimOpBoundaryProcess::GetWrappedProcess() const
{
  return static_cast<G4OpBoundaryProcess*>(pRegProcess);
}
