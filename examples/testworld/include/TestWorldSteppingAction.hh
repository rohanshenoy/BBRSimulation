#ifndef TestWorldSteppingAction_hh
#define TestWorldSteppingAction_hh

#include "G4OpBoundaryProcess.hh"
#include "G4UserSteppingAction.hh"

class TestWorldRunAction;
class BBSimOpBoundaryProcess;

// Fills the G4Analysis "crossings" ntuple per optical-photon boundary crossing
// and the "abspoints" ntuple on photon termination. The ntuples and the
// categorical encoders are owned by TestWorldRunAction (non-owning pointer here).
class TestWorldSteppingAction : public G4UserSteppingAction {
public:
  explicit TestWorldSteppingAction(TestWorldRunAction* runAction);
  ~TestWorldSteppingAction() override = default;
  void UserSteppingAction(const G4Step*) override;

private:
  TestWorldRunAction*            fRunAction;        // non-owning; lifetime managed by ActionInit
  BBSimOpBoundaryProcess*  fWrapper        = nullptr;
  G4OpBoundaryProcess*     fBoundary       = nullptr;
  G4int                    fCurrentRunID   = -1;
  G4int                    fCurrentTrackID = -1;
  G4int                    fCurrentEventID = -1;
  G4int                    fNReflect       = 0;
};

#endif
