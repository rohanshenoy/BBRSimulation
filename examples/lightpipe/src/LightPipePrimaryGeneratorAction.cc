#include "LightPipePrimaryGeneratorAction.hh"

void LightPipePrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fSource.Generate(event);
}
