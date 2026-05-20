#include "BBRTestActionInit.hh"
#include "BBRRunAction.hh"
#include "BBRTestPGA.hh"
#include "BBRTestSteppingAction.hh"

// MT master: RunAction only (no stepping, no primary generation on master).
void BBRTestActionInit::BuildForMaster() const
{
  SetUserAction(new BBRRunAction());
}

// Workers (and sole thread in serial mode): all actions.
void BBRTestActionInit::Build() const
{
  auto* runAction = new BBRRunAction();
  SetUserAction(runAction);
  SetUserAction(new BBRTestPGA());
  SetUserAction(new BBRTestSteppingAction(runAction));
}
