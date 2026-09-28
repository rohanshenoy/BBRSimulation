#include "TestWorldActionInitialization.hh"
#include "TestWorldRunAction.hh"
#include "TestWorldPrimaryGeneratorAction.hh"
#include "TestWorldSteppingAction.hh"
#include "BBRConfigManager.hh"

// MT master: RunAction only (no stepping, no primary generation on master).
void TestWorldActionInitialization::BuildForMaster() const
{
  SetUserAction(new TestWorldRunAction());
}

// Workers (and sole thread in serial mode): all actions.
void TestWorldActionInitialization::Build() const
{
  BBRConfigManager::Instance();   // create this worker's config clone + messenger early

  auto* runAction = new TestWorldRunAction();
  SetUserAction(runAction);
  SetUserAction(new TestWorldPrimaryGeneratorAction());
  SetUserAction(new TestWorldSteppingAction(runAction));
}
