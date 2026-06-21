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

  fSetTCmd = new G4UIcmdWithADouble("/bbr/thermal/setT", this);
  fSetTCmd->SetGuidance("Planck emitter temperature [K].");
  fSetTCmd->SetParameterName("T", false);
  fSetTCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunModeCmd = new G4UIcmdWithABool("/bbr/gun/mode", this);
  fGunModeCmd->SetGuidance("true = fixed gun, false = Planck emitter.");
  fGunModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunPosXCmd = new G4UIcmdWithADouble("/bbr/gun/posX", this);
  fGunPosXCmd->SetGuidance("Gun X position [mm].");
  fGunPosXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunPosYCmd = new G4UIcmdWithADouble("/bbr/gun/posY", this);
  fGunPosYCmd->SetGuidance("Gun Y position [mm].");
  fGunPosYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunPosZCmd = new G4UIcmdWithADouble("/bbr/gun/posZ", this);
  fGunPosZCmd->SetGuidance("Gun Z position [mm].");
  fGunPosZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunDirXCmd = new G4UIcmdWithADouble("/bbr/gun/dirX", this);
  fGunDirXCmd->SetGuidance("Gun direction X component.");
  fGunDirXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunDirYCmd = new G4UIcmdWithADouble("/bbr/gun/dirY", this);
  fGunDirYCmd->SetGuidance("Gun direction Y component.");
  fGunDirYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunDirZCmd = new G4UIcmdWithADouble("/bbr/gun/dirZ", this);
  fGunDirZCmd->SetGuidance("Gun direction Z component.");
  fGunDirZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGunECmd = new G4UIcmdWithADouble("/bbr/gun/energy_eV", this);
  fGunECmd->SetGuidance("Photon energy [eV] (500 GHz = 2.07e-3 eV).");
  fGunECmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCuMatCmd = new G4UIcmdWithAString("/bbr/det/setCuMaterial", this);
  fCuMatCmd->SetGuidance("Named Cu alias: OFHC_Cu (RRR=100) | OF_Cu (RRR=3) | "
                         "HP_Cu (RRR=6). Sets RRR and resets stage T to 4 K.");
  fCuMatCmd->SetParameterName("alias", true);
  fCuMatCmd->SetDefaultValue("OFHC_Cu");
  fCuMatCmd->SetToBeBroadcasted(false);
  fCuMatCmd->AvailableForStates(G4State_PreInit);

  fCuRRRCmd = new G4UIcmdWithAnInteger("/bbr/det/setCuRRR", this);
  fCuRRRCmd->SetGuidance("Cu Residual Resistance Ratio (integer >= 1). "
                         "sigma_DC = RRR * 5.96e7 S/m. Before /run/initialize.");
  fCuRRRCmd->SetParameterName("RRR", false);
  fCuRRRCmd->SetToBeBroadcasted(false);
  fCuRRRCmd->AvailableForStates(G4State_PreInit);

  fCuStageTCmd = new G4UIcmdWithADoubleAndUnit("/bbr/det/setCuStageT", this);
  fCuStageTCmd->SetGuidance("Temperature stage [K] for the Cu reflectance table "
                            "(default 4 K). Before /run/initialize.");
  fCuStageTCmd->SetParameterName("T", false);
  fCuStageTCmd->SetDefaultUnit("kelvin");
  fCuStageTCmd->SetUnitCandidates("kelvin");
  fCuStageTCmd->SetToBeBroadcasted(false);
  fCuStageTCmd->AvailableForStates(G4State_PreInit);
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
  if      (cmd == fPrintCmd)   { BBRConfigManager::Print(G4cout); }
  else if (cmd == fSetTCmd)    { BBRConfigManager::SetThermalT_K(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunModeCmd) { BBRConfigManager::SetGunMode(G4UIcmdWithABool::GetNewBoolValue(value)); }
  else if (cmd == fGunPosXCmd) { BBRConfigManager::SetGunPosX_mm(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunPosYCmd) { BBRConfigManager::SetGunPosY_mm(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunPosZCmd) { BBRConfigManager::SetGunPosZ_mm(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunDirXCmd) { BBRConfigManager::SetGunDirX(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunDirYCmd) { BBRConfigManager::SetGunDirY(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunDirZCmd) { BBRConfigManager::SetGunDirZ(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fGunECmd)    { BBRConfigManager::SetGunEnergy_eV(G4UIcmdWithADouble::GetNewDoubleValue(value)); }
  else if (cmd == fCuMatCmd)   { BBRConfigManager::SetCuMaterial(value); }
  else if (cmd == fCuRRRCmd)   { BBRConfigManager::SetCuRRR(G4UIcmdWithAnInteger::GetNewIntValue(value)); }
  else if (cmd == fCuStageTCmd){ BBRConfigManager::SetCuStageT_K(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(value)); }
}
