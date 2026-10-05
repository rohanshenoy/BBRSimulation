#include "BBRConfigManager.hh"
#include "BBRConfigMessenger.hh"

#include "G4Exception.hh"
#include "G4Threading.hh"
#include <cmath>
#include <cstdlib>
#include "G4ios.hh"
#include <ostream>

#ifndef BBRSIM_DEFAULT_DATADIR
#error "BBRSIM_DEFAULT_DATADIR must be defined by the build (library/CMakeLists.txt)"
#endif

BBRConfigManager* BBRConfigManager::Instance() {
  static const BBRConfigManager* masterInstance = nullptr;
  static G4ThreadLocal BBRConfigManager* theInstance = nullptr;

  if (!theInstance) {
    if (!G4Threading::IsWorkerThread()) {        // master or sequential
      theInstance = new BBRConfigManager;
      masterInstance = theInstance;
    } else {                                     // workers clone from master
      if (!masterInstance) G4Exception("BBRConfigManager::Instance", "BBR020", FatalException, "a worker thread asked for the config before the master created it; call BBRConfigManager::Instance() in main() first");
      theInstance = new BBRConfigManager(*masterInstance);
    }
  }
  return theInstance;
}

BBRConfigManager::BBRConfigManager()
  : fThermalT_K(4.0),
    fEmitterCenter_mm(-50.0, 0.0, 0.0),
    fEmitterSize_mm(1.0, 20.0, 20.0),
    fGunMode(false),
    fGunPosX_mm(-20.0), fGunPosY_mm(0.0), fGunPosZ_mm(0.0),
    fGunDirX(1.0), fGunDirY(0.0), fGunDirZ(0.0),
    fGunEnergy_eV(2.07e-3),
    fGunPol(0.0, 0.0, 0.0),
    fDataDir(std::getenv("BBRSIMDATA") ? std::getenv("BBRSIMDATA") : BBRSIM_DEFAULT_DATADIR),
    fCuRRR(100), fCuStageT_K(4.0),
    fMessenger(new BBRConfigMessenger(this)) {}

BBRConfigManager::BBRConfigManager(const BBRConfigManager& master)
  : fThermalT_K(master.fThermalT_K),
    fEmitterCenter_mm(master.fEmitterCenter_mm),
    fEmitterSize_mm(master.fEmitterSize_mm),
    fGunMode(master.fGunMode),
    fGunPosX_mm(master.fGunPosX_mm), fGunPosY_mm(master.fGunPosY_mm),
    fGunPosZ_mm(master.fGunPosZ_mm),
    fGunDirX(master.fGunDirX), fGunDirY(master.fGunDirY),
    fGunDirZ(master.fGunDirZ),
    fGunEnergy_eV(master.fGunEnergy_eV),
    fGunPol(master.fGunPol),
    fDataDir(master.fDataDir),
    fCuRRR(master.fCuRRR), fCuStageT_K(master.fCuStageT_K),
    fMessenger(new BBRConfigMessenger(this)) {}

BBRConfigManager::~BBRConfigManager() { delete fMessenger; fMessenger = nullptr; }

G4bool BBRConfigManager::SetThermalT_K(G4double v) {
  if (!std::isfinite(v) || v <= 0.) {
    G4cerr << "[BBR] thermal/setT: temperature must be > 0 K, got " << v << G4endl;
    return false;
  }
  Instance()->fThermalT_K = v;
  return true;
}

G4bool BBRConfigManager::SetEmitterSize_mm(const G4ThreeVector& v) {
  if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z()) ||
      v.x() <= 0. || v.y() <= 0. || v.z() <= 0.) {
    G4cerr << "[BBR] thermal/emitterSize: all extents must be > 0, got "
           << v << " mm" << G4endl;
    return false;
  }
  Instance()->fEmitterSize_mm = v;
  return true;
}

G4bool BBRConfigManager::SetGunEnergy_eV(G4double v) {
  if (!std::isfinite(v) || v <= 0.) {
    G4cerr << "[BBR] gun/energy_eV: energy must be > 0 eV, got " << v << G4endl;
    return false;
  }
  Instance()->fGunEnergy_eV = v;
  return true;
}

G4bool BBRConfigManager::SetCuRRR(G4int rrr) {
  if (rrr < 1) {
    G4cerr << "[BBR] det/setCuRRR: RRR must be >= 1, got " << rrr << G4endl;
    return false;
  }
  Instance()->fCuRRR = rrr;
  return true;
}

G4bool BBRConfigManager::SetCuStageT_K(G4double T_K) {
  if (!std::isfinite(T_K) || T_K <= 0.) {
    G4cerr << "[BBR] det/setCuStageT: temperature must be > 0 K, got " << T_K << G4endl;
    return false;
  }
  Instance()->fCuStageT_K = T_K;
  return true;
}

G4bool BBRConfigManager::SetCuMaterial(const G4String& alias) {
  auto* m = Instance();
  if      (alias == "OFHC_Cu") { m->fCuRRR = 100; m->fCuStageT_K = 4.0; }
  else if (alias == "OF_Cu")   { m->fCuRRR =   3; m->fCuStageT_K = 4.0; }
  else if (alias == "HP_Cu")   { m->fCuRRR =   6; m->fCuStageT_K = 4.0; }
  else {
    G4cerr << "[BBR] det/setCuMaterial: unknown alias '" << alias
           << "'.  Valid: OFHC_Cu, OF_Cu, HP_Cu." << G4endl;
    return false;
  }
  return true;
}

void BBRConfigManager::printConfig(std::ostream& os) const {
  os << "=== BBRConfigManager settings ===\n"
     << "  /bbr/thermal/setT     " << fThermalT_K   << " K\n"
     << "  /bbr/thermal/emitterCenter " << fEmitterCenter_mm << " mm\n"
     << "  /bbr/thermal/emitterSize   " << fEmitterSize_mm   << " mm\n"
     << "  /bbr/gun/mode         " << (fGunMode ? "true" : "false") << "\n"
     << "  /bbr/gun/pos[XYZ]     " << fGunPosX_mm << " " << fGunPosY_mm << " "
                                   << fGunPosZ_mm << " mm\n"
     << "  /bbr/gun/dir[XYZ]     " << fGunDirX << " " << fGunDirY << " "
                                   << fGunDirZ << "\n"
     << "  /bbr/gun/energy_eV    " << fGunEnergy_eV << " eV\n"
     << "  /bbr/gun/pol          " << fGunPol
     << (fGunPol.mag2() > 0. ? "" : "  (random)") << "\n"
     << "  /bbr/dataDir          " << fDataDir << "   # BBRSIMDATA\n"
     << "  /bbr/det/setCuRRR     " << fCuRRR << "\n"
     << "  /bbr/det/setCuStageT  " << fCuStageT_K  << " K\n"
     << "=================================" << G4endl;
}
