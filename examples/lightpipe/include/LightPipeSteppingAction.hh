#ifndef LightPipeSteppingAction_hh
#define LightPipeSteppingAction_hh

#include "G4OpBoundaryProcess.hh"
#include "G4UserSteppingAction.hh"

class LightPipeRunAction;
class BBSimOpBoundaryProcess;

// Fills the G4Analysis "crossings" ntuple per optical-photon boundary crossing
// and the "abspoints" ntuple on photon termination. The ntuples and the
// categorical encoders are owned by LightPipeRunAction (non-owning pointer here).
class LightPipeSteppingAction : public G4UserSteppingAction {
public:
  explicit LightPipeSteppingAction(LightPipeRunAction* runAction);
  ~LightPipeSteppingAction() override = default;
  void UserSteppingAction(const G4Step*) override;

private:
  LightPipeRunAction*            fRunAction;        // non-owning; lifetime managed by ActionInit
  BBSimOpBoundaryProcess*  fWrapper        = nullptr;
  G4OpBoundaryProcess*     fBoundary       = nullptr;
  G4int                    fCurrentRunID   = -1;
  G4int                    fCurrentTrackID = -1;
  G4int                    fCurrentEventID = -1;
  G4int                    fNReflect       = 0;
};

#endif
