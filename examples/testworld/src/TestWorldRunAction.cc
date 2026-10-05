#include "TestWorldRunAction.hh"
#include "BBRApplicationInfo.hh"
TestWorldRunAction::TestWorldRunAction()
  : fAnalysis(BBRSIM_APP_BUILD_VERSION, BBRSIM_APP_SOURCE_FINGERPRINT) {}
void TestWorldRunAction::BeginOfRunAction(const G4Run* run) { fAnalysis.BeginRun(run, IsMaster()); }
void TestWorldRunAction::EndOfRunAction(const G4Run* run) { fAnalysis.EndRun(run, IsMaster()); }
