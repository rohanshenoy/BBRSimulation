#ifndef TestWorldRunAction_hh
#define TestWorldRunAction_hh
#include "G4UserRunAction.hh"
#include "BBRAnalysis.hh"

// Each example owns its actions; the ordinary helper owns output bookkeeping.
class TestWorldRunAction : public G4UserRunAction {
 public:
  TestWorldRunAction();
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;
  BBRAnalysis& Analysis() { return fAnalysis; }
 private:
  BBRAnalysis fAnalysis;
};
#endif
