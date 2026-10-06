#include "BBSimOpBoundaryProcess.hh"

#include "BBRCrackLibrary.hh"
#include "G4AffineTransform.hh"
#include "G4Exception.hh"
#include "G4GeometryTolerance.hh"
#include "G4NavigationHistory.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SafetyHelper.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"

#include <cmath>

BBSimOpBoundaryProcess::BBSimOpBoundaryProcess(const G4String& name)
  : G4WrapperProcess(name)
{}

void BBSimOpBoundaryProcess::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
  G4WrapperProcess::BuildPhysicsTable(particle);
  BBRCrackLibrary::Instance().ValidatePlacedCracks();
}

// PostStepDoIt — intercept steps that ENTER a vacuum_wg crack volume (HFSS
// diffraction) or a material carrying a REFLECTIVITY table (tabulated
// reflectance); fall through to the stock process otherwise. Geometries with
// neither are handled identically to stock G4OpBoundaryProcess.
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
  fLastHFSSFreqGHz = -1.;

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

  // Dataset for this crack volume at this photon's frequency (lazy-loaded).
  // Units: in Geant4 internals E is in MeV and h_Planck in MeV*ns, so
  // E/h_Planck is a frequency in 1/ns; 1e9*hertz is exactly 1/ns, so the
  // quotient is in GHz.
  const G4String volName   = touch->GetVolume()->GetName();
  const G4String datasetId = BBRCrackLibrary::DatasetIdOf(volName);
  const G4double nu_GHz =
      aTrack.GetKineticEnergy() / CLHEP::h_Planck / (1e9 * CLHEP::hertz);
  const BBRHFSSData& hfss =
      BBRCrackLibrary::Instance().Lookup(datasetId, nu_GHz, fLastHFSSFreqGHz);

  // --- incoming angles in crack-local frame (folded into HFSS quarter-symmetry) ---
  // HFSS global frame in crack-local axes: X = -phi_hat (gap), Y = +theta_hat
  // (long), Z = +normal_hat (propagation). HFSS's spherical incidence angles are
  // those of the ARRIVAL direction r of the plane wave, k = -r:
  //   IWaveTheta = polar angle of -k from Z   = acos(-k.normal_hat)
  //   IWavePhi   = azimuth of -k from X to Y  = atan2(-k.theta_hat, k.phi_hat)
  // A photon entering along normal_hat has IWaveTheta = 180 deg (normal entry;
  // crack1's raw T is 1.055 there and ~0 at 0 deg). The exit coordinate system
  // (x_e, y_e, z_e) = (Z, Y, -X) is then exactly (normal_hat, theta_hat, phi_hat),
  // so the CSV exit positions and far-field angles need no transform.
  G4double cosVal = std::min(1., std::max(-1., -khat.dot(normal_hat)));
  G4double iwaveTheta_deg = std::acos(cosVal) * (180. / CLHEP::pi);
  // IWavePhi: azimuth of -k in the HFSS X-Y plane (X = -phi_hat, Y = +theta_hat). Transverse
  // components below 1e-12 are snapped to +0 first: for a k exactly in the
  // x-z plane (k_y == 0) atan2(-0, -x) is -180° but atan2(+0, -x) is +180°,
  // and which one the expression yields depends on how the optimiser orders
  // the dot product. The two results mirror theta_hat (sy) differently; both
  // are physically equivalent by the plate's y-symmetry, but the choice must
  // be deterministic so the Python mirror (bbrsim.hfss) reproduces it.
  G4double kt = -khat.dot(theta_hat);
  G4double kp =  khat.dot(phi_hat);
  if (std::abs(kt) < 1e-12) kt = 0.;
  if (std::abs(kp) < 1e-12) kp = 0.;
  G4double iwavePhi_raw   = std::atan2(kt, kp) * (180. / CLHEP::pi);

  // The HFSS sweep covers IWavePhi ∈ [0°, 90°] only; the parallel-plate
  // geometry is mirror-symmetric about both transverse axes. Fold the azimuth
  // into that wedge by mirroring the transverse axes, and carry the mirror
  // signs so everything sampled from the folded dataset (outgoing direction,
  // polarization, exit position) is mapped back to the true frame. Without
  // the un-fold, oblique photons get mirror-image outgoing distributions.
  //   sy: phi → |phi|        (mirrors Y, i.e. flips theta_hat)
  //   sx: |phi| → 180−|phi|  (mirrors X, i.e. flips phi_hat)
  const G4double sy = (iwavePhi_raw < 0.) ? -1. : 1.;
  G4double iwavePhi_deg = std::abs(iwavePhi_raw);
  const G4double sx = (iwavePhi_deg > 90.) ? -1. : 1.;
  if (iwavePhi_deg > 90.) iwavePhi_deg = 180. - iwavePhi_deg;

  // Folded-frame transverse axes expressed in world coordinates. Building all
  // basis vectors from these applies the mirror to inputs (polarization
  // decomposition) and un-applies it to sampled outputs in one stroke.
  const G4ThreeVector theta_f = sy * theta_hat;
  const G4ThreeVector phi_f   = sx * phi_hat;

  // --- decompose incoming polarization onto the HFSS incoming basis ---
  // In the folded HFSS frame (X, Y, Z) = (-phi_f, theta_f, normal_hat) the
  // spherical basis at the arrival direction r(T, P) is
  //   e_theta = cosT cosP X + cosT sinP Y - sinT Z,   e_phi = -sinP X + cosP Y.
  // The two vectors below are exactly -e_theta and -e_phi. The common sign
  // cancels in T, in the cross term and in |E|^2, and flips pol_out only, which
  // is the same state. Do not flip the +sin(th) term on its own: it is the
  // -e_theta component, not a sign error.
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

  // --- transmittance decision (Wang eq. 58 in eq. 53, with the cross term) ---
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

    // The HFSS model is non-local: it samples the state at the far face of
    // the crack. Keep the track just inside that same crack volume and update
    // Geant4's navigator before proposing the displacement. The following
    // transportation step will then cross crack -> World normally and update
    // the track touchable/material state. Placing the track directly outside
    // would leave the ordinary G4ParticleChange touchable inconsistent.
    const G4double surfaceTolerance =
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    const G4double inset = 10. * surfaceTolerance;
    if (halfLen <= inset) {
      G4Exception("BBSimOpBoundaryProcess::HandleDiffractionBoundary",
                  "BBR004", FatalException,
                  "Crack half-length is too small for the navigation inset.");
    }

    const G4ThreeVector pos_inside = pos_out - inset * normal_hat;
    const G4ThreeVector local_inside = xf.TransformPoint(pos_inside);
    if (touch->GetSolid()->Inside(local_inside) != kInside) {
      G4Exception("BBSimOpBoundaryProcess::HandleDiffractionBoundary",
                  "BBR005", FatalException,
                  "HFSS exit sample is not inside the crack volume.");
    }

    auto* safetyHelper =
        G4TransportationManager::GetTransportationManager()->GetSafetyHelper();
    safetyHelper->Locate(pos_inside, dir_out);

    fParticleChange.ProposeMomentumDirection(dir_out);
    fParticleChange.ProposePolarization(pol_out);
    fParticleChange.ProposePosition(pos_inside);
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

  // Surface normal of the boundary just crossed, taken from the navigator
  // exactly as stock G4OpBoundaryProcess does (theGlobalNormal). It points out
  // of the exited volume into mat2; flip it to oppose k (back into the vacuum).
  //
  // Do NOT derive it from the entered solid's SurfaceNormal() in general: when
  // the entered volume is the MOTHER of the exited one (a photon inside a
  // vacuum_wg crack striking the crack's side wall and entering the Cu slab),
  // the hit point is interior to the slab box and G4Box::SurfaceNormal
  // silently returns the nearest slab face (x) instead of the wall normal (z),
  // which reflects the photon about the wrong axis and sends it through solid
  // copper. The entered solid's normal is kept only as a defensive fallback
  // (no state tried so far reaches it: the navigator reported a valid normal in
  // every case, including the crack-wall one) for
  // the rare case where the navigator cannot provide one; that fallback is
  // correct whenever the entered solid's surface IS the boundary (the usual
  // vacuum -> metal hit), so it is reported as a warning, not a fatal error.
  const G4ThreeVector hitPos = aStep.GetPostStepPoint()->GetPosition();
  G4bool validNormal = false;
  const G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
  auto iNav = G4TransportationManager::GetTransportationManager()
                  ->GetActiveNavigatorsIterator();
  G4ThreeVector nhat = (iNav[hNavId])->GetGlobalExitNormal(hitPos, &validNormal);
  if (!validNormal) {
    const G4VTouchable*      postTouch = aStep.GetPostStepPoint()->GetTouchable();
    const G4AffineTransform& postXF   = postTouch->GetHistory()->GetTopTransform();
    const G4ThreeVector      posLocal = postXF.TransformPoint(hitPos);
    nhat = postXF.InverseTransformAxis(postTouch->GetSolid()->SurfaceNormal(posLocal));
    G4ExceptionDescription ed;
    ed << "Navigator returned no valid exit normal at " << hitPos / mm << " mm ("
       << aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() << " -> "
       << aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()
       << "); falling back to the entered solid's SurfaceNormal(), which is "
          "wrong if the photon is leaving a daughter volume into this one.";
    G4Exception("BBSimOpBoundaryProcess::HandleReflectanceBoundary", "BBR006",
                JustWarning, ed);
  }
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
