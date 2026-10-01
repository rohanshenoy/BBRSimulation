#include "TestWorldSteppingAction.hh"
#include "TestWorldRunAction.hh"
TestWorldSteppingAction::TestWorldSteppingAction(TestWorldRunAction* runAction)
  : fRecorder(&runAction->Analysis()) {}
void TestWorldSteppingAction::UserSteppingAction(const G4Step* step) { fRecorder.Record(step); }
