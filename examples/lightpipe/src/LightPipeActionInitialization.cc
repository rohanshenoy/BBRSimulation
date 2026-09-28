#include "LightPipeActionInitialization.hh"
#include "LightPipeRunAction.hh"
#include "LightPipePrimaryGeneratorAction.hh"
#include "LightPipeSteppingAction.hh"
#include "BBRConfigManager.hh"

// MT master: RunAction only (no stepping, no primary generation on master).
void LightPipeActionInitialization::BuildForMaster() const
{
  SetUserAction(new LightPipeRunAction());
}

// Workers (and sole thread in serial mode): all actions.
void LightPipeActionInitialization::Build() const
{
  BBRConfigManager::Instance();   // create this worker's config clone + messenger early

  auto* runAction = new LightPipeRunAction();
  SetUserAction(runAction);
  SetUserAction(new LightPipePrimaryGeneratorAction());
  SetUserAction(new LightPipeSteppingAction(runAction));
}
