#ifndef LightPipeActionInitialization_hh
#define LightPipeActionInitialization_hh

#include "G4VUserActionInitialization.hh"

// Wires all BBR test user actions.
// BuildForMaster: registers a LightPipeRunAction for the master thread (MT only).
// Build: registers LightPipeRunAction + LightPipePrimaryGeneratorAction + LightPipeSteppingAction.
// The run manager is multithreaded (G4RunManagerFactory); in a sequential
// build only Build() is invoked.
class LightPipeActionInitialization : public G4VUserActionInitialization {
public:
  LightPipeActionInitialization()           = default;
  ~LightPipeActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build()          const override;
};

#endif
