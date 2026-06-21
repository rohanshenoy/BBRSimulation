#ifndef BBRConfigMessenger_hh
#define BBRConfigMessenger_hh

#include "G4UImessenger.hh"

class BBRConfigManager;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcommand;

// Hand-written messenger owning every /bbr/... command and forwarding to
// BBRConfigManager. One messenger is constructed per thread-local manager
// (see BBRConfigManager copy ctor), so broadcast commands reach worker copies.
class BBRConfigMessenger : public G4UImessenger {
 public:
  explicit BBRConfigMessenger(BBRConfigManager* mgr);
  ~BBRConfigMessenger() override;

  void SetNewValue(G4UIcommand* cmd, G4String value) override;

 private:
  BBRConfigManager* fManager;

  G4UIcmdWithoutParameter*   fPrintCmd   = nullptr;

  // thermal + gun (wired in Task 2)
  G4UIcmdWithADouble*        fSetTCmd    = nullptr;
  G4UIcmdWithABool*          fGunModeCmd = nullptr;
  G4UIcmdWithADouble*        fGunPosXCmd = nullptr;
  G4UIcmdWithADouble*        fGunPosYCmd = nullptr;
  G4UIcmdWithADouble*        fGunPosZCmd = nullptr;
  G4UIcmdWithADouble*        fGunDirXCmd = nullptr;
  G4UIcmdWithADouble*        fGunDirYCmd = nullptr;
  G4UIcmdWithADouble*        fGunDirZCmd = nullptr;
  G4UIcmdWithADouble*        fGunECmd    = nullptr;

  // det (wired in Task 3)
  G4UIcmdWithAString*        fCuMatCmd   = nullptr;
  G4UIcmdWithAnInteger*      fCuRRRCmd   = nullptr;
  G4UIcmdWithADoubleAndUnit* fCuStageTCmd = nullptr;
};

#endif
