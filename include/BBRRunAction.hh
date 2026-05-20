#ifndef BBRRunAction_hh
#define BBRRunAction_hh

#include "G4UserRunAction.hh"
#include <fstream>

// Owns test_output.csv for the duration of a run.
// BBRTestSteppingAction writes rows via GetOutputStream().
// Note: single-threaded (serial RunManager) only. MT would require
// either per-thread filenames or a G4Run::Merge accumulation design.
class BBRRunAction : public G4UserRunAction {
public:
  BBRRunAction();
  ~BBRRunAction() override;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run)   override;

  std::ofstream& GetOutputStream() { return fOut; }

private:
  std::ofstream fOut;
};

#endif
