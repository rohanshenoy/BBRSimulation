#ifndef BBRTestSteppingAction_hh
#define BBRTestSteppingAction_hh

#include "G4OpBoundaryProcess.hh"
#include "G4UserSteppingAction.hh"
#include <fstream>

class BBRTestSteppingAction : public G4UserSteppingAction {
public:
  BBRTestSteppingAction();
  ~BBRTestSteppingAction() override;
  void UserSteppingAction(const G4Step*) override;

private:
  std::ofstream        fOut;
  G4OpBoundaryProcess* fBoundary       = nullptr;
  G4int                fCurrentTrackID = -1;
  G4int                fNReflect       = 0;
};
#endif
