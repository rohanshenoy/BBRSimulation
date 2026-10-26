#include "TestWorldMessenger.hh"
#include "TestWorldDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

#include <sstream>

TestWorldMessenger::TestWorldMessenger(TestWorldDetectorConstruction* det) : fDet(det)
{
  fDir = new G4UIdirectory("/bbr/testworld/");
  fDir->SetGuidance("Test-world geometry options (before /run/initialize).");
  fRoundGapCmd = new G4UIcmdWithABool("/bbr/testworld/roundGap", this);
  using DC = TestWorldDetectorConstruction;
  std::ostringstream line1, line2;
  line1 << "Also place the straight round gap RoundGap_r50um: a " << DC::kRoundGapLength / mm
        << " mm Cu plate at";
  line2 << "z = " << DC::kRoundGapPlateZ / mm << " mm with a " << DC::kRoundGapRadius / um
        << " um-radius vacuum_wg hole along x (HFSS radius " << DC::kRoundGapHfssRadius / um
        << " um).";
  fRoundGapCmd->SetGuidance(line1.str().c_str());
  fRoundGapCmd->SetGuidance(line2.str().c_str());
  fRoundGapCmd->SetParameterName("enable", false);
  fRoundGapCmd->AvailableForStates(G4State_PreInit);
  fRoundGapCmd->SetToBeBroadcasted(false);

  fCracksCmd = new G4UIcmdWithABool("/bbr/testworld/cracks", this);
  fCracksCmd->SetGuidance("Place the Cu plates CuPlateCrack1 and CuPlateCrack2 with their vacuum_wg");
  fCracksCmd->SetGuidance("cracks (default true); false leaves them out. The round gap is unaffected.");
  fCracksCmd->SetParameterName("enable", false);
  fCracksCmd->AvailableForStates(G4State_PreInit);
  fCracksCmd->SetToBeBroadcasted(false);
}

TestWorldMessenger::~TestWorldMessenger()
{
  delete fCracksCmd;
  delete fRoundGapCmd;
  delete fDir;
}

void TestWorldMessenger::SetNewValue(G4UIcommand* cmd, G4String value)
{
  if (cmd == fRoundGapCmd) fDet->SetRoundGap(G4UIcmdWithABool::GetNewBoolValue(value));
  if (cmd == fCracksCmd) fDet->SetCracks(G4UIcmdWithABool::GetNewBoolValue(value));
}
