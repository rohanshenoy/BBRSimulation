#ifndef TestWorldActionInitialization_hh
#define TestWorldActionInitialization_hh

#include "G4VUserActionInitialization.hh"

// Wires all BBR test user actions.
// BuildForMaster: registers a TestWorldRunAction for the master thread (MT only).
// Build: registers TestWorldRunAction + TestWorldPrimaryGeneratorAction + TestWorldSteppingAction.
// The run manager is multithreaded (G4RunManagerFactory); in a sequential
// build only Build() is invoked.
class TestWorldActionInitialization : public G4VUserActionInitialization {
public:
  TestWorldActionInitialization()           = default;
  ~TestWorldActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build()          const override;
};

#endif
