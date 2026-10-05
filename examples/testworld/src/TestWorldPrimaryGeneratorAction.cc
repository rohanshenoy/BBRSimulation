#include "TestWorldPrimaryGeneratorAction.hh"

void TestWorldPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fSource.Generate(event);
}
