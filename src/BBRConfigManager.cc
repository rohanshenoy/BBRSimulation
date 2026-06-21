#include "BBRConfigManager.hh"
#include "BBRConfigMessenger.hh"

#include "G4Threading.hh"
#include "G4ios.hh"
#include <ostream>

BBRConfigManager* BBRConfigManager::Instance() {
  static const BBRConfigManager* masterInstance = nullptr;
  static G4ThreadLocal BBRConfigManager* theInstance = nullptr;

  if (!theInstance) {
    if (!G4Threading::IsWorkerThread()) {        // master or sequential
      theInstance = new BBRConfigManager;
      masterInstance = theInstance;
    } else {                                     // workers clone from master
      theInstance = new BBRConfigManager(*masterInstance);
    }
  }
  return theInstance;
}

BBRConfigManager::BBRConfigManager()
  : fThermalT_K(4.0),
    fGunMode(false),
    fGunPosX_mm(-20.0), fGunPosY_mm(0.0), fGunPosZ_mm(0.0),
    fGunDirX(1.0), fGunDirY(0.0), fGunDirZ(0.0),
    fGunEnergy_eV(2.07e-3),
    fCuRRR(100), fCuStageT_K(4.0),
    fMessenger(new BBRConfigMessenger(this)) {}

BBRConfigManager::BBRConfigManager(const BBRConfigManager& master)
  : fThermalT_K(master.fThermalT_K),
    fGunMode(master.fGunMode),
    fGunPosX_mm(master.fGunPosX_mm), fGunPosY_mm(master.fGunPosY_mm),
    fGunPosZ_mm(master.fGunPosZ_mm),
    fGunDirX(master.fGunDirX), fGunDirY(master.fGunDirY),
    fGunDirZ(master.fGunDirZ),
    fGunEnergy_eV(master.fGunEnergy_eV),
    fCuRRR(master.fCuRRR), fCuStageT_K(master.fCuStageT_K),
    fMessenger(new BBRConfigMessenger(this)) {}

BBRConfigManager::~BBRConfigManager() { delete fMessenger; fMessenger = nullptr; }

void BBRConfigManager::SetThermalT_K(G4double v) {
  if (v <= 0.) {
    G4cerr << "[BBR] thermal/setT: temperature must be > 0 K, got " << v << G4endl;
    return;
  }
  Instance()->fThermalT_K = v;
}

void BBRConfigManager::SetCuRRR(G4int rrr) {
  if (rrr < 1) {
    G4cerr << "[BBR] det/setCuRRR: RRR must be >= 1, got " << rrr << G4endl;
    return;
  }
  Instance()->fCuRRR = rrr;
}

void BBRConfigManager::SetCuStageT_K(G4double T_K) {
  if (T_K <= 0.) {
    G4cerr << "[BBR] det/setCuStageT: temperature must be > 0 K, got " << T_K << G4endl;
    return;
  }
  Instance()->fCuStageT_K = T_K;
}

void BBRConfigManager::SetCuMaterial(const G4String& alias) {
  auto* m = Instance();
  if      (alias == "OFHC_Cu") { m->fCuRRR = 100; m->fCuStageT_K = 4.0; }
  else if (alias == "OF_Cu")   { m->fCuRRR =   3; m->fCuStageT_K = 4.0; }
  else if (alias == "HP_Cu")   { m->fCuRRR =   6; m->fCuStageT_K = 4.0; }
  else {
    G4cerr << "[BBR] det/setCuMaterial: unknown alias '" << alias
           << "'.  Valid: OFHC_Cu, OF_Cu, HP_Cu." << G4endl;
  }
}

void BBRConfigManager::printConfig(std::ostream& os) const {
  os << "=== BBRConfigManager settings ===\n"
     << "  /bbr/thermal/setT     " << fThermalT_K   << " K\n"
     << "  /bbr/gun/mode         " << (fGunMode ? "true" : "false") << "\n"
     << "  /bbr/gun/pos[XYZ]     " << fGunPosX_mm << " " << fGunPosY_mm << " "
                                   << fGunPosZ_mm << " mm\n"
     << "  /bbr/gun/dir[XYZ]     " << fGunDirX << " " << fGunDirY << " "
                                   << fGunDirZ << "\n"
     << "  /bbr/gun/energy_eV    " << fGunEnergy_eV << " eV\n"
     << "  /bbr/det/setCuRRR     " << fCuRRR << "\n"
     << "  /bbr/det/setCuStageT  " << fCuStageT_K  << " K\n"
     << "=================================\n";
}
