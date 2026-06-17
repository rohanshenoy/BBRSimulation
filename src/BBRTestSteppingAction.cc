#include "BBRTestSteppingAction.hh"
#include "BBRRunAction.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4AnalysisManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <cmath>

namespace {
G4String StatusStr(G4OpBoundaryProcessStatus s) {
  switch (s) {
    case FresnelRefraction:       return "FresnelRefraction";
    case FresnelReflection:       return "FresnelReflection";
    case TotalInternalReflection: return "TIR";
    case LambertianReflection:    return "LambertianReflection";
    case LobeReflection:          return "LobeReflection";
    case SpikeReflection:         return "SpikeReflection";
    case BackScattering:          return "BackScattering";
    case Absorption:              return "Absorption";
    case Detection:               return "Detection";
    case NotAtBoundary:           return "NotAtBoundary";
    case SameMaterial:            return "SameMaterial";
    case StepTooSmall:            return "StepTooSmall";
    case NoRINDEX:                return "NoRINDEX";
    default:                      return "Other";
  }
}

// Resolve the boundary outcome for the current step. If the BBR wrapper handled
// the boundary itself, the wrapped stock process never ran (its GetStatus()
// would be stale), so prefer the wrapper's own record; fall back to the stock
// status, then to "unknown".
G4String ResolveStatus(BBSimOpBoundaryProcess* wrapper,
                       G4OpBoundaryProcess* boundary) {
  if (wrapper) {
    G4String s = wrapper->GetLastBBRStatusString();
    if (!s.empty()) return s;
  }
  if (boundary) return StatusStr(boundary->GetStatus());
  return "unknown";
}
} // namespace

BBRTestSteppingAction::BBRTestSteppingAction(BBRRunAction* runAction)
  : G4UserSteppingAction()
  , fRunAction(runAction)
{}

void BBRTestSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhoton()) return;

  const G4int runId   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  const G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // Per-track reflection counter — reset when run, event, or track changes.
  // Done for every optical-photon step (not only boundary steps) so a photon
  // that terminates without ever crossing a boundary (e.g. future bulk
  // absorption in a dielectric) still reports a per-track-correct count in the
  // termination block below. The increment stays gated to real crossings.
  if (runId != fCurrentRunID || eventId != fCurrentEventID ||
      track->GetTrackID() != fCurrentTrackID) {
    fCurrentRunID   = runId;
    fCurrentEventID = eventId;
    fCurrentTrackID = track->GetTrackID();
    fNReflect = 0;
  }

  // Record a termination point (one row per killed optical photon). n_reflect
  // here is the number of boundary crossings BEFORE the terminating contact.
  if (track->GetTrackStatus() == fStopAndKill) {
    const G4StepPoint* end = step->GetPostStepPoint();
    const G4ThreeVector ep = end->GetPosition();
    const G4ThreeVector ed = end->GetMomentumDirection();
    const G4String tvol =
        end->GetPhysicalVolume() ? end->GetPhysicalVolume()->GetName() : "none";
    const G4String tstat = ResolveStatus(fWrapper, fBoundary);
    auto* am = G4AnalysisManager::Instance();
    const auto& a = fRunAction->fAbs;
    const G4int aid = fRunAction->fAbsPointsId;
    am->FillNtupleIColumn(aid, a.run_id, runId);
    am->FillNtupleIColumn(aid, a.event_id, eventId);
    am->FillNtupleDColumn(aid, a.x, ep.x() / mm);
    am->FillNtupleDColumn(aid, a.y, ep.y() / mm);
    am->FillNtupleDColumn(aid, a.z, ep.z() / mm);
    am->FillNtupleDColumn(aid, a.energy, end->GetKineticEnergy() / eV);
    am->FillNtupleDColumn(aid, a.px, ed.x());
    am->FillNtupleDColumn(aid, a.py, ed.y());
    am->FillNtupleDColumn(aid, a.pz, ed.z());
    am->FillNtupleIColumn(aid, a.n_reflect, fNReflect);
    am->FillNtupleIColumn(aid, a.term_vol, fRunAction->EncodeVolume(tvol));
    am->FillNtupleIColumn(aid, a.term_status, fRunAction->EncodeStatus(tstat));
    am->AddNtupleRow(aid);
  }

  if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary)  return;

  // Skip tolerance-scale re-steps at the same surface (StepTooSmall in the
  // boundary process): no physics happens on them and they would inflate
  // n_reflect / duplicate rows.
  const G4double kCarTolerance =
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  if (step->GetStepLength() <= kCarTolerance / 2.) return;

  // Lazy-init: find BBSimOpBoundaryProcess wrapper on first boundary step
  if (!fWrapper) {
    G4ProcessVector* pv = track->GetDefinition()
                               ->GetProcessManager()->GetProcessList();
    for (std::size_t i = 0; i < pv->size(); ++i) {
      if (auto* w = dynamic_cast<BBSimOpBoundaryProcess*>((*pv)[i])) {
        fWrapper  = w;
        fBoundary = w->GetWrappedProcess();
        break;
      }
    }
  }

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  G4ThreeVector      pPre = pre->GetMomentumDirection();

  ++fNReflect;  // count this boundary crossing (counter reset per-track above)

  const G4ThreeVector pos   = post->GetPosition();
  const G4ThreeVector pPost = post->GetMomentumDirection();

  const G4String volPre  = pre->GetPhysicalVolume()
                             ? pre->GetPhysicalVolume()->GetName()  : "none";
  const G4String matPre  = pre->GetMaterial()
                             ? pre->GetMaterial()->GetName()         : "none";
  const G4String volPost = post->GetPhysicalVolume()
                             ? post->GetPhysicalVolume()->GetName() : "none";
  const G4String matPost = post->GetMaterial()
                             ? post->GetMaterial()->GetName()        : "none";

  const G4String status = ResolveStatus(fWrapper, fBoundary);

  const G4double theta_in = std::acos(std::abs(pPre.x())) * 180. / CLHEP::pi;
  const G4double phi_in = std::atan2(pPre.y(), pPre.z()) * 180. / CLHEP::pi;

  auto* am = G4AnalysisManager::Instance();
  const auto& c = fRunAction->fCross;
  const G4int id = fRunAction->fCrossingsId;
  am->FillNtupleIColumn(id, c.run_id, runId);
  am->FillNtupleIColumn(id, c.event_id, eventId);
  am->FillNtupleDColumn(id, c.x, pos.x() / mm);
  am->FillNtupleDColumn(id, c.y, pos.y() / mm);
  am->FillNtupleDColumn(id, c.z, pos.z() / mm);
  am->FillNtupleDColumn(id, c.energy, pre->GetKineticEnergy() / eV);
  am->FillNtupleDColumn(id, c.px_pre, pPre.x());
  am->FillNtupleDColumn(id, c.py_pre, pPre.y());
  am->FillNtupleDColumn(id, c.pz_pre, pPre.z());
  am->FillNtupleDColumn(id, c.px_post, pPost.x());
  am->FillNtupleDColumn(id, c.py_post, pPost.y());
  am->FillNtupleDColumn(id, c.pz_post, pPost.z());
  am->FillNtupleDColumn(id, c.theta_in, theta_in);
  am->FillNtupleDColumn(id, c.phi_in, phi_in);
  am->FillNtupleIColumn(id, c.vol_pre, fRunAction->EncodeVolume(volPre));
  am->FillNtupleIColumn(id, c.mat_pre, fRunAction->EncodeMaterial(matPre));
  am->FillNtupleIColumn(id, c.vol_post, fRunAction->EncodeVolume(volPost));
  am->FillNtupleIColumn(id, c.mat_post, fRunAction->EncodeMaterial(matPost));
  am->FillNtupleIColumn(id, c.status, fRunAction->EncodeStatus(status));
  am->FillNtupleIColumn(id, c.event_type,
                        BBRRunAction::EventTypeForStatus(status));
  am->FillNtupleIColumn(id, c.n_reflect, fNReflect);
  am->AddNtupleRow(id);
}
