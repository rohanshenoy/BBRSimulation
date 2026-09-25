#ifndef BBSimPhysics_hh
#define BBSimPhysics_hh

#include "G4VPhysicsConstructor.hh"

class G4ProcessManager;

// Physics constructor that finds the G4OpBoundaryProcess registered by
// G4OpticalPhysics, removes it, and replaces it with BBSimOpBoundaryProcess.
// Must be registered AFTER G4OpticalPhysics so the boundary process exists.
class BBSimPhysics : public G4VPhysicsConstructor
{
 public:
  explicit BBSimPhysics(G4int verbose = 0);
  ~BBSimPhysics() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;

 private:
  void WrapOpBoundaryProcess(G4ProcessManager* procMan) const;
};

#endif
