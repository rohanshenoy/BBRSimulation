#ifndef TestWorldMessenger_hh
#define TestWorldMessenger_hh

#include "G4UImessenger.hh"

class G4UIcmdWithABool;
class G4UIdirectory;
class TestWorldDetectorConstruction;

// /bbr/testworld/ commands, owned by the test world (not the library):
//   /bbr/testworld/roundGap <bool>   also place the straight round gap
//                                    RoundGap_r50um (default false)
// PreInit only and not broadcast: the geometry is built on the master.
class TestWorldMessenger : public G4UImessenger {
 public:
  explicit TestWorldMessenger(TestWorldDetectorConstruction* det);
  ~TestWorldMessenger() override;
  void SetNewValue(G4UIcommand* cmd, G4String value) override;

 private:
  TestWorldDetectorConstruction* fDet;
  G4UIdirectory* fDir;
  G4UIcmdWithABool* fRoundGapCmd;
};

#endif
