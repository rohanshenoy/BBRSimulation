#include "BBRTestActionInit.hh"
#include "BBRRunAction.hh"
#include "BBRTestPGA.hh"
#include "BBRTestSteppingAction.hh"
#include "BBRConfigManager.hh"

// MT master: RunAction only (no stepping, no primary generation on master).
void BBRTestActionInit::BuildForMaster() const
{
  SetUserAction(new BBRRunAction());
}

// Workers (and sole thread in serial mode): all actions.
void BBRTestActionInit::Build() const
{
  BBRConfigManager::Instance();   // create this worker's config clone + messenger early

  auto* runAction = new BBRRunAction();
  SetUserAction(runAction);
  SetUserAction(new BBRTestPGA());
  SetUserAction(new BBRTestSteppingAction(runAction));
}
