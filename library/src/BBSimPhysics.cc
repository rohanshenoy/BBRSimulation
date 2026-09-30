#include "BBSimPhysics.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include <vector>

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

// Find the G4OpBoundaryProcess in the optical photon process manager, remove
// it, wrap it with BBSimOpBoundaryProcess and put the wrapper in its place.
// Follows the pattern of CDMSRDecayPhysics::WrapRDMProcess() (SuperSim).
// BBR014 (fatal): a wrapper is already registered (BBSimPhysics registered
// twice) or there is no stock process to wrap (registered before
// G4OpticalPhysics, which would then add a second, unwrapped boundary process).
void BBSimPhysics::WrapOpBoundaryProcess(G4ProcessManager* procMan) const
{
  if(!procMan) return;
  G4ProcessVector* procs = procMan->GetProcessList();
  std::size_t iBoundary = procs->size();
  for(std::size_t i = 0; i < procs->size(); ++i) {
    if(dynamic_cast<BBSimOpBoundaryProcess*>((*procs)[i]))
      G4Exception("BBSimPhysics::WrapOpBoundaryProcess", "BBR014", FatalException,
                  "G4OpBoundaryProcess is already wrapped: BBSimPhysics was registered twice.");
    if(iBoundary == procs->size() && dynamic_cast<G4OpBoundaryProcess*>((*procs)[i])) iBoundary = i;
  }
  if(iBoundary == procs->size())
    G4Exception("BBSimPhysics::WrapOpBoundaryProcess", "BBR014", FatalException,
                "No G4OpBoundaryProcess to wrap: register BBSimPhysics after G4OpticalPhysics.");

  // Put the wrapper exactly where the stock process was. Every discrete process
  // draws a random number when it resets its interaction length, so a wrapper
  // appended after OpWLS/OpWLS2 would change which process gets which number, and
  // a run with WLS active would no longer be bit-identical to stock. All optical
  // processes share the default ordering parameter, so an ordering value cannot
  // express the position (AddDiscreteProcess(p, ord) appends after equal
  // orderings): take the processes that follow the stock one off the list, insert
  // the wrapper, and put them back in their original order with their orderings.
  struct Entry { G4VProcess* proc; G4int ordAtRest, ordAlong, ordPost; };
  auto entryOf = [procMan](G4VProcess* p) {
    return Entry{p, procMan->GetProcessOrdering(p, idxAtRest),
                 procMan->GetProcessOrdering(p, idxAlongStep),
                 procMan->GetProcessOrdering(p, idxPostStep)};
  };
  const Entry stock = entryOf((*procs)[iBoundary]);
  std::vector<Entry> after;
  for(std::size_t i = iBoundary + 1; i < procs->size(); ++i) after.push_back(entryOf((*procs)[i]));
  for(const auto& e : after) procMan->RemoveProcess(e.proc);
  procMan->RemoveProcess(stock.proc);

  auto* wrapper = new BBSimOpBoundaryProcess();
  wrapper->RegisterProcess(stock.proc);
  procMan->AddProcess(wrapper, stock.ordAtRest, stock.ordAlong, stock.ordPost);
  for(const auto& e : after) procMan->AddProcess(e.proc, e.ordAtRest, e.ordAlong, e.ordPost);
}
