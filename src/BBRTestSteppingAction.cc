#include "BBRTestSteppingAction.hh"
#include "BBRRunAction.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
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

  // Lazy-init: find BBSimOpBoundaryProcess wrapper on first boundary step
  if (!fBoundary) {
    G4ProcessVector* pv = track->GetDefinition()
                               ->GetProcessManager()->GetProcessList();
    for (std::size_t i = 0; i < pv->size(); ++i) {
      if (auto* w = dynamic_cast<BBSimOpBoundaryProcess*>((*pv)[i])) {
        fBoundary = dynamic_cast<G4OpBoundaryProcess*>(w->GetWrappedProcess());
        break;
      }
    }
  }

  G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  G4ThreeVector      pPre = pre->GetMomentumDirection();

  // Skip StepTooSmall steps — pre-step momentum direction is undefined
  if (std::abs(pPre.x()) > 1.1 || std::abs(pPre.y()) > 1.1 || std::abs(pPre.z()) > 1.1) return;

  // Per-track reflection counter — reset when event or track changes
  if (eventId != fCurrentEventID || track->GetTrackID() != fCurrentTrackID) {
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
  const G4String status  = fBoundary ? StatusStr(fBoundary->GetStatus()) : "unknown";

  const G4double theta_in = std::acos(std::abs(pPre.x())) * 180. / CLHEP::pi;
  const G4double phi_in   = std::atan2(pPre.y(), pPre.z()) * 180. / CLHEP::pi;

  std::ofstream& out = fRunAction->GetOutputStream();
  out << eventId << ","
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
