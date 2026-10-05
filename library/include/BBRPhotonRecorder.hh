#ifndef BBRPhotonRecorder_hh
#define BBRPhotonRecorder_hh

#include "G4OpBoundaryProcess.hh"
#include "G4Step.hh"

class BBRAnalysis;
class BBSimOpBoundaryProcess;

// Fills the G4Analysis "crossings" ntuple per optical-photon boundary crossing
// and the "abspoints" ntuple on photon termination. The ntuples and the
// categorical encoders are owned by BBRAnalysis (non-owning pointer here).
class BBRPhotonRecorder {
public:
  explicit BBRPhotonRecorder(BBRAnalysis* runAction);
  ~BBRPhotonRecorder() = default;
  void Record(const G4Step*);

private:
  BBRAnalysis*            fRunAction;        // non-owning; lifetime managed by ActionInit
  BBSimOpBoundaryProcess*  fWrapper        = nullptr;
  G4OpBoundaryProcess*     fBoundary       = nullptr;
  G4int                    fCurrentRunID   = -1;
  G4int                    fCurrentTrackID = -1;
  G4int                    fCurrentEventID = -1;
  G4int                    fNReflect       = 0;
  G4int                    fReflections    = 0;
};

#endif
