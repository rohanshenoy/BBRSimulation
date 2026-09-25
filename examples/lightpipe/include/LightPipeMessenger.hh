#ifndef LightPipeMessenger_hh
#define LightPipeMessenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"

class LightPipeDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

// Owns the /bbr/lightpipe/ geometry commands (G4CMP per-construction messenger
// pattern; NOT folded into BBRConfigManager). All PreInit + non-broadcast:
// geometry is built on the master in Construct().
class LightPipeMessenger : public G4UImessenger {
public:
  explicit LightPipeMessenger(LightPipeDetectorConstruction* det);
  ~LightPipeMessenger() override;

  void SetNewValue(G4UIcommand* cmd, G4String val) override;

private:
  LightPipeDetectorConstruction* fDet;
  G4UIdirectory*             fDir;
  G4UIcmdWithADoubleAndUnit* fBoreCmd;
  G4UIcmdWithADoubleAndUnit* fLengthCmd;
  G4UIcmdWithADoubleAndUnit* fWallCmd;
  G4UIcmdWithAString*        fWallMatCmd;
  G4UIcmdWithAString*        fModeCmd;
  G4UIcmdWithAString*        fStlPathCmd;
};

#endif
