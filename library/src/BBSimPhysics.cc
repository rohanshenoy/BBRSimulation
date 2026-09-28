#include "BBSimPhysics.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

BBSimPhysics::BBSimPhysics(G4int verbose)
  : G4VPhysicsConstructor("BBSimPhysics")
{
  SetVerboseLevel(verbose);
}

void BBSimPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhoton();
}

void BBSimPhysics::ConstructProcess()
{
  if(verboseLevel > 1)
    G4cout << "BBSimPhysics::ConstructProcess()" << G4endl;

  WrapOpBoundaryProcess(
    G4OpticalPhoton::OpticalPhoton()->GetProcessManager());
}

// Find G4OpBoundaryProcess in the optical photon process manager, remove it,
// wrap it with BBSimOpBoundaryProcess, and re-add the wrapper.
// Follows the pattern of CDMSRDecayPhysics::WrapRDMProcess() (SuperSim).
void BBSimPhysics::WrapOpBoundaryProcess(G4ProcessManager* procMan) const
{
  if(!procMan) return;

  G4VProcess* opBoundary = nullptr;
  G4ProcessVector* procs = procMan->GetProcessList();
  if(procs)
  {
    for(size_t i = 0; i < procs->size() && !opBoundary; ++i)
    {
      if(dynamic_cast<G4OpBoundaryProcess*>((*procs)[i]))
        opBoundary = (*procs)[i];
    }
  }

  if(opBoundary)
  {
    procMan->RemoveProcess(opBoundary);
  }
  else
  {
    if(verboseLevel > 0)
      G4cerr << "BBSimPhysics WARNING: G4OpBoundaryProcess not registered; "
             << "creating a new one." << G4endl;
    opBoundary = new G4OpBoundaryProcess();
  }

  auto* wrapper = new BBSimOpBoundaryProcess();
  wrapper->RegisterProcess(opBoundary);
  procMan->AddDiscreteProcess(wrapper);
}
