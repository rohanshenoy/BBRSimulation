#include "LightPipeRunAction.hh"
#include "BBRApplicationInfo.hh"
LightPipeRunAction::LightPipeRunAction()
  : fAnalysis(BBRSIM_APP_BUILD_VERSION, BBRSIM_APP_SOURCE_FINGERPRINT) {}
void LightPipeRunAction::BeginOfRunAction(const G4Run* run) { fAnalysis.BeginRun(run, IsMaster()); }
void LightPipeRunAction::EndOfRunAction(const G4Run* run) { fAnalysis.EndRun(run, IsMaster()); }
