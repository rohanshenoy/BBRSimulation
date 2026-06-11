#ifndef BBRTestSteppingAction_hh
#define BBRTestSteppingAction_hh

#include "G4OpBoundaryProcess.hh"
#include "G4UserSteppingAction.hh"
#include <fstream>

class BBRRunAction;
class BBSimOpBoundaryProcess;

// Records every optical photon boundary crossing to test_output.csv.
// The output stream is owned by BBRRunAction; this class holds a non-owning pointer.
class BBRTestSteppingAction : public G4UserSteppingAction {
public:
  explicit BBRTestSteppingAction(BBRRunAction* runAction);
  ~BBRTestSteppingAction() override = default;
  void UserSteppingAction(const G4Step*) override;

private:
  BBRRunAction*            fRunAction;        // non-owning; lifetime managed by ActionInit
  BBSimOpBoundaryProcess*  fWrapper        = nullptr;
  G4OpBoundaryProcess*     fBoundary       = nullptr;
  G4int                    fCurrentRunID   = -1;
  G4int                    fCurrentTrackID = -1;
  G4int                    fCurrentEventID = -1;
  G4int                    fNReflect       = 0;
};

#endif
