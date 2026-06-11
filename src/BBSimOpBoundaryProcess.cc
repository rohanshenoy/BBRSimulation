#include "BBSimOpBoundaryProcess.hh"

#include "BBRCrackLibrary.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
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
//
// The wrapped G4OpBoundaryProcess is Forced, so this runs on EVERY step of
// every optical photon. Mirror the stock process's own entry guards before
// intercepting: only act on genuine geometry-boundary steps longer than the
// surface tolerance. Tolerance-scale re-steps at the same surface must NOT
// roll the absorption dice a second time — the stock process classifies them
// as StepTooSmall and does nothing.
G4VParticleChange* BBSimOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                                        const G4Step& aStep)
{
  fLastBBRStatus = kBBRNone;

  const G4double kCarTolerance =
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary &&
      aStep.GetStepLength() > kCarTolerance / 2.) {
    const G4Material* mat2 = aStep.GetPostStepPoint()->GetMaterial();
    if (mat2 && mat2->GetName() == "vacuum_wg")
      return HandleDiffractionBoundary(aTrack, aStep);
    if (mat2) {
      G4MaterialPropertiesTable* mpt = mat2->GetMaterialPropertiesTable();
      if (mpt && mpt->GetProperty("REFLECTIVITY"))
        return HandleReflectanceBoundary(aTrack, aStep);
    }
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

  // The HFSS sweep covers IWavePhi ∈ [0°, 90°] only; the parallel-plate
  // geometry is mirror-symmetric about both transverse axes. Fold the azimuth
  // into that wedge by mirroring the transverse axes, and carry the mirror
  // signs so everything sampled from the folded dataset (outgoing direction,
  // polarization, exit position) is mapped back to the true frame. Without
  // the un-fold, oblique photons get mirror-image outgoing distributions.
  //   sy: phi → |phi|        (mirrors ŷ_i, i.e. flips theta_hat)
  //   sx: |phi| → 180−|phi|  (mirrors x̂_i, i.e. flips phi_hat)
  const G4double sy = (iwavePhi_raw < 0.) ? -1. : 1.;
  G4double iwavePhi_deg = std::abs(iwavePhi_raw);
  const G4double sx = (iwavePhi_deg > 90.) ? -1. : 1.;
  if (iwavePhi_deg > 90.) iwavePhi_deg = 180. - iwavePhi_deg;

  // Folded-frame transverse axes expressed in world coordinates. Building all
  // basis vectors from these applies the mirror to inputs (polarization
  // decomposition) and un-applies it to sampled outputs in one stroke.
  const G4ThreeVector theta_f = sy * theta_hat;
  const G4ThreeVector phi_f   = sx * phi_hat;

  // --- decompose incoming polarization onto HFSS incoming spherical basis ---
  // Standard spherical basis at (IWaveTheta, IWavePhi) in the folded HFSS
  // frame (ẑ_i=normal_hat, x̂_i=phi_f, ŷ_i=-theta_f):
  //   ê_θ = -sin(T)*ẑ_i + cos(T)*cos(P)*x̂_i - cos(T)*sin(P)*ŷ_i
  //   ê_φ =            -sin(P)*x̂_i           - cos(P)*ŷ_i
  G4double th = iwaveTheta_deg * (CLHEP::pi / 180.);
  G4double ph = iwavePhi_deg   * (CLHEP::pi / 180.);
  G4ThreeVector eTheta_in = +std::sin(th) * normal_hat
                            - std::cos(th) * std::sin(ph) * theta_f
                            + std::cos(th) * std::cos(ph) * phi_f;
  G4ThreeVector ePhi_in   = -std::cos(ph) * theta_f
                            - std::sin(ph) * phi_f;

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

  if (transmitted) {
    fLastBBRStatus = kBBRDiffractionTransmit;

    // Transmit: sample outgoing direction + polarization (Wang eqs. 56-57).
    // The folded axes (theta_f, phi_f) un-mirror the sampled output.
    G4ThreeVector pol_out;
    G4ThreeVector dir_out = hfss.SampleOutgoingDirection(
        E_theta, E_phi, iwavePhi_deg, iwaveTheta_deg,
        phi_f, theta_f, normal_hat, pol_out);

    // Sample exit-face position (Wang eq. 55).
    // Exit face center: project the volume's local origin onto the exit face
    // along the local outward normal. The HFSS waveguide CSV stores exit
    // positions as ABSOLUTE cross-section coordinates centred on the
    // waveguide axis, so they must be added to the face centre. (Adding them
    // to the entry-point projection — the previous behaviour — relocated
    // off-axis photons outside the crack, into solid Cu.)
    G4ThreeVector localNormal = xf.TransformAxis(normal_hat);
    G4double      halfLen     = touch->GetSolid()->DistanceToOut(
                                    G4ThreeVector(0., 0., 0.), localNormal);
    G4ThreeVector exitCenter  = xf.InverseTransformPoint(halfLen * localNormal);

    // HFSS waveguide CSV: Y=long dim → theta axis, Z=gap/b → phi axis.
    G4ThreeVector pos_out = hfss.SampleExitPosition(
        E_theta, E_phi, iwavePhi_deg, iwaveTheta_deg,
        exitCenter, theta_f, phi_f);

    fParticleChange.ProposeMomentumDirection(dir_out);
    fParticleChange.ProposePolarization(pol_out);
    fParticleChange.ProposePosition(pos_out);
    fParticleChange.ProposeTrackStatus(fAlive);
  } else {
    fLastBBRStatus = kBBRDiffractionReflect;

    // Reflect: specular flip of momentum and polarization about the crack normal.
    G4ThreeVector dir_ref = (khat - 2.*khat.dot(normal_hat)*normal_hat).unit();
    G4ThreeVector pol_ref = phat  - 2.*phat.dot(normal_hat)*normal_hat;
    if (pol_ref.mag() > 1e-30) pol_ref = pol_ref.unit();
    else                        pol_ref = phi_hat;

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
    fLastBBRStatus = kBBRAbsorb;
    ++fNReflectanceAbsorb;
    fParticleChange.ProposeLocalEnergyDeposit(E);
    fParticleChange.ProposeTrackStatus(fStopAndKill);
  } else {
    fLastBBRStatus = kBBRReflect;
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

G4String BBSimOpBoundaryProcess::GetLastBBRStatusString() const
{
  switch (fLastBBRStatus) {
    case kBBRDiffractionTransmit: return "BBRDiffractionTransmit";
    case kBBRDiffractionReflect:  return "BBRDiffractionReflect";
    case kBBRReflect:             return "BBRReflect";
    case kBBRAbsorb:              return "BBRAbsorb";
    case kBBRNone: default:       return "";
  }
}

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
