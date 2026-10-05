#ifndef LightPipeRunAction_hh
#define LightPipeRunAction_hh
#include "G4UserRunAction.hh"
#include "BBRAnalysis.hh"

// Each example owns its actions; the ordinary helper owns output bookkeeping.
class LightPipeRunAction : public G4UserRunAction {
 public:
  LightPipeRunAction();
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;
  BBRAnalysis& Analysis() { return fAnalysis; }
 private:
  BBRAnalysis fAnalysis;
};
#endif
