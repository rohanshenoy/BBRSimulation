#ifndef BBRTestActionInit_hh
#define BBRTestActionInit_hh

#include "G4VUserActionInitialization.hh"

// Wires all BBR test user actions.
// BuildForMaster: registers a BBRRunAction for the master thread (MT only).
// Build: registers BBRRunAction + BBRTestPGA + BBRTestSteppingAction.
// In serial mode (current default) only Build() is invoked.
class BBRTestActionInit : public G4VUserActionInitialization {
public:
  BBRTestActionInit()           = default;
  ~BBRTestActionInit() override = default;

  void BuildForMaster() const override;
  void Build()          const override;
};

#endif
