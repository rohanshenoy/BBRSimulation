#include "BBRTestSteppingAction.hh"
#include "BBRRunAction.hh"
#include "BBSimOpBoundaryProcess.hh"

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
#include <iomanip>

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
} // namespace

BBRTestSteppingAction::BBRTestSteppingAction(BBRRunAction* runAction)
  : G4UserSteppingAction()
  , fRunAction(runAction)
{}

void BBRTestSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhoton()) return;
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

  G4int runId   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  G4ThreeVector      pPre = pre->GetMomentumDirection();

  // Per-track reflection counter — reset when run, event, or track changes
  if (runId != fCurrentRunID || eventId != fCurrentEventID ||
      track->GetTrackID() != fCurrentTrackID) {
    fCurrentRunID   = runId;
    fCurrentEventID = eventId;
    fCurrentTrackID = track->GetTrackID();
    fNReflect = 0;
  }
  ++fNReflect;

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

  // Status: if the wrapper handled this boundary itself the wrapped stock
  // process never ran, so its GetStatus() would be stale — prefer the
  // wrapper's own record of what it did.
  G4String status = "unknown";
  if (fWrapper) {
    status = fWrapper->GetLastBBRStatusString();
    if (status.empty() && fBoundary) status = StatusStr(fBoundary->GetStatus());
  } else if (fBoundary) {
    status = StatusStr(fBoundary->GetStatus());
  }

  const G4double theta_in = std::acos(std::abs(pPre.x())) * 180. / CLHEP::pi;
  const G4double phi_in   = std::atan2(pPre.y(), pPre.z()) * 180. / CLHEP::pi;

  std::ofstream& out = fRunAction->GetOutputStream();
  out << runId << "," << eventId << ","
      << std::fixed << std::setprecision(6)
      << pos.x()/mm << "," << pos.y()/mm << "," << pos.z()/mm << ","
      << pre->GetKineticEnergy()/eV << ","
      << std::setprecision(10)
      << pPre.x()  << "," << pPre.y()  << "," << pPre.z()  << ","
      << pPost.x() << "," << pPost.y() << "," << pPost.z() << ","
      << std::setprecision(6)
      << theta_in << "," << phi_in << ","
      << volPre << "," << matPre << ","
      << volPost << "," << matPost << ","
      << status << "," << fNReflect << "\n";
}
