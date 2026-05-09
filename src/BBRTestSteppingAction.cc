#include "BBRTestSteppingAction.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
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

BBRTestSteppingAction::BBRTestSteppingAction()
  : G4UserSteppingAction()
{
  fOut.open("test_output.csv");
  fOut << "event_id,x_mm,y_mm,z_mm,energy_eV,"
          "px_pre,py_pre,pz_pre,px_post,py_post,pz_post,"
          "vol_pre,mat_pre,vol_post,mat_post,status,n_reflect\n";
}

BBRTestSteppingAction::~BBRTestSteppingAction()
{
  if (fOut.is_open()) fOut.close();
}

void BBRTestSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhoton()) return;
  if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary)  return;

  // Lazy-init: find G4OpBoundaryProcess via BBSimOpBoundaryProcess wrapper
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

  // Per-track reflection counter
  if (track->GetTrackID() != fCurrentTrackID) {
    fCurrentTrackID = track->GetTrackID();
    fNReflect = 0;
  }
  ++fNReflect;

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  G4ThreeVector      pos  = post->GetPosition();
  G4ThreeVector      pPre = pre->GetMomentumDirection();
  G4ThreeVector      pPost= post->GetMomentumDirection();

  G4String volPre  = pre->GetPhysicalVolume()
                       ? pre->GetPhysicalVolume()->GetName()  : "none";
  G4String matPre  = pre->GetMaterial()
                       ? pre->GetMaterial()->GetName()         : "none";
  G4String volPost = post->GetPhysicalVolume()
                       ? post->GetPhysicalVolume()->GetName() : "none";
  G4String matPost = post->GetMaterial()
                       ? post->GetMaterial()->GetName()        : "none";
  G4String status  = fBoundary ? StatusStr(fBoundary->GetStatus()) : "unknown";

  G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  fOut << eventId << ","
       << std::fixed << std::setprecision(6)
       << pos.x()/mm << "," << pos.y()/mm << "," << pos.z()/mm << ","
       << pre->GetKineticEnergy()/eV << ","
       << pPre.x()  << "," << pPre.y()  << "," << pPre.z()  << ","
       << pPost.x() << "," << pPost.y() << "," << pPost.z() << ","
       << volPre << "," << matPre << ","
       << volPost << "," << matPost << ","
       << status << "," << fNReflect << "\n";
}
