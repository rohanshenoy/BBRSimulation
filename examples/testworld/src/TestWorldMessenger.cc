#include "TestWorldMessenger.hh"
#include "TestWorldDetectorConstruction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

TestWorldMessenger::TestWorldMessenger(TestWorldDetectorConstruction* det) : fDet(det)
{
  fDir = new G4UIdirectory("/bbr/testworld/");
  fDir->SetGuidance("Test-world geometry options (before /run/initialize).");
  fRoundGapCmd = new G4UIcmdWithABool("/bbr/testworld/roundGap", this);
  fRoundGapCmd->SetGuidance("Also place the straight round gap RoundGap_r50um: a 0.4 mm Cu plate at");
  fRoundGapCmd->SetGuidance("z = -80 mm with a 51 um-radius vacuum_wg hole along x (HFSS radius 50 um).");
  fRoundGapCmd->SetParameterName("enable", false);
  fRoundGapCmd->AvailableForStates(G4State_PreInit);
  fRoundGapCmd->SetToBeBroadcasted(false);
}

TestWorldMessenger::~TestWorldMessenger()
{
  delete fRoundGapCmd;
  delete fDir;
}

void TestWorldMessenger::SetNewValue(G4UIcommand* cmd, G4String value)
{
  if (cmd == fRoundGapCmd) fDet->SetRoundGap(G4UIcmdWithABool::GetNewBoolValue(value));
}
