#include "LightPipeSteppingAction.hh"
#include "LightPipeRunAction.hh"
LightPipeSteppingAction::LightPipeSteppingAction(LightPipeRunAction* runAction)
  : fRecorder(&runAction->Analysis()) {}
void LightPipeSteppingAction::UserSteppingAction(const G4Step* step) { fRecorder.Record(step); }
