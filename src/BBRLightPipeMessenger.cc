#include "BBRLightPipeMessenger.hh"
#include "BBRLightPipeDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

BBRLightPipeMessenger::BBRLightPipeMessenger(
    BBRLightPipeDetectorConstruction* det)
  : fDet(det)
{
  fDir = new G4UIdirectory("/bbr/lightpipe/");
  fDir->SetGuidance("Light-pipe geometry parameters (PreInit only).");

  fBoreCmd = new G4UIcmdWithADoubleAndUnit("/bbr/lightpipe/bore", this);
  fBoreCmd->SetGuidance("Inner bore radius (aperture).");
  fBoreCmd->SetParameterName("bore", false);
  fBoreCmd->SetUnitCategory("Length");
  fBoreCmd->AvailableForStates(G4State_PreInit);
  fBoreCmd->SetToBeBroadcasted(false);

  fLengthCmd = new G4UIcmdWithADoubleAndUnit("/bbr/lightpipe/length", this);
  fLengthCmd->SetGuidance("Tube length along +x.");
  fLengthCmd->SetParameterName("length", false);
  fLengthCmd->SetUnitCategory("Length");
  fLengthCmd->AvailableForStates(G4State_PreInit);
  fLengthCmd->SetToBeBroadcasted(false);

  fWallCmd = new G4UIcmdWithADoubleAndUnit("/bbr/lightpipe/wallThickness", this);
  fWallCmd->SetGuidance("Wall thickness.");
  fWallCmd->SetParameterName("wallThickness", false);
  fWallCmd->SetUnitCategory("Length");
  fWallCmd->AvailableForStates(G4State_PreInit);
  fWallCmd->SetToBeBroadcasted(false);

  fWallMatCmd = new G4UIcmdWithAString("/bbr/lightpipe/wallMaterial", this);
  fWallMatCmd->SetGuidance("Wall optical material: Cu | reflector.");
  fWallMatCmd->SetParameterName("wallMaterial", false);
  fWallMatCmd->SetCandidates("Cu reflector");
  fWallMatCmd->AvailableForStates(G4State_PreInit);
  fWallMatCmd->SetToBeBroadcasted(false);
}

BBRLightPipeMessenger::~BBRLightPipeMessenger()
{
  delete fBoreCmd;
  delete fLengthCmd;
  delete fWallCmd;
  delete fWallMatCmd;
  delete fDir;
}

void BBRLightPipeMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
  if (cmd == fBoreCmd)
    fDet->SetBore(fBoreCmd->GetNewDoubleValue(val));
  else if (cmd == fLengthCmd)
    fDet->SetLength(fLengthCmd->GetNewDoubleValue(val));
  else if (cmd == fWallCmd)
    fDet->SetWallThickness(fWallCmd->GetNewDoubleValue(val));
  else if (cmd == fWallMatCmd)
    fDet->SetWallMaterial(val);
}
