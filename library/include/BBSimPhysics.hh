#ifndef BBSimPhysics_hh
#define BBSimPhysics_hh

#include "G4VPhysicsConstructor.hh"

class G4ProcessManager;

// Physics constructor that finds the G4OpBoundaryProcess registered by
// G4OpticalPhysics, removes it, and replaces it with BBSimOpBoundaryProcess at
// the same place in the optical photon's process list (the processes that
// followed it keep their order), so a geometry without vacuum_wg volumes or
// material REFLECTIVITY runs bit-identically to stock, also with WLS active.
// Must be registered AFTER G4OpticalPhysics, and once: otherwise
// WrapOpBoundaryProcess raises the fatal BBR014.
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
