#ifndef LightPipeSteppingAction_hh
#define LightPipeSteppingAction_hh
#include "G4UserSteppingAction.hh"
#include "BBRPhotonRecorder.hh"
class LightPipeRunAction;
class LightPipeSteppingAction : public G4UserSteppingAction {
 public:
  explicit LightPipeSteppingAction(LightPipeRunAction* runAction);
  void UserSteppingAction(const G4Step* step) override;
 private:
  BBRPhotonRecorder fRecorder;
};
#endif
