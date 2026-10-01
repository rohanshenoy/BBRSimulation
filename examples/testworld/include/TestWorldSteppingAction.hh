#ifndef TestWorldSteppingAction_hh
#define TestWorldSteppingAction_hh
#include "G4UserSteppingAction.hh"
#include "BBRPhotonRecorder.hh"
class TestWorldRunAction;
class TestWorldSteppingAction : public G4UserSteppingAction {
 public:
  explicit TestWorldSteppingAction(TestWorldRunAction* runAction);
  void UserSteppingAction(const G4Step* step) override;
 private:
  BBRPhotonRecorder fRecorder;
};
#endif
