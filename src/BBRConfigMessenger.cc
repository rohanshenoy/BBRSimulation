#include "BBRConfigMessenger.hh"
#include "BBRConfigManager.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"

BBRConfigMessenger::BBRConfigMessenger(BBRConfigManager* mgr)
  : G4UImessenger(), fManager(mgr) {
  fPrintCmd = new G4UIcmdWithoutParameter("/bbr/config/print", this);
  fPrintCmd->SetGuidance("Print all BBRConfigManager settings.");
  fPrintCmd->SetToBeBroadcasted(false);
  // (thermal/gun/det commands registered in later tasks)
}

BBRConfigMessenger::~BBRConfigMessenger() {
  delete fPrintCmd;
  delete fSetTCmd;
  delete fGunModeCmd;
  delete fGunPosXCmd; delete fGunPosYCmd; delete fGunPosZCmd;
  delete fGunDirXCmd; delete fGunDirYCmd; delete fGunDirZCmd;
  delete fGunECmd;
  delete fCuMatCmd; delete fCuRRRCmd; delete fCuStageTCmd;
}

void BBRConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == fPrintCmd) { BBRConfigManager::Print(G4cout); return; }
  (void)value;
}
