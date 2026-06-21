#ifndef BBRConfigManager_hh
#define BBRConfigManager_hh

// Centralized, thread-local run configuration for BBRsim.
//
// Pattern mirrors G4CMP's G4CMPConfigManager: Instance() returns a per-thread
// object. The master/sequential thread constructs the real instance (compiled
// defaults); each worker thread copy-constructs from the master and owns its
// own BBRConfigMessenger, so broadcast commands reach every worker's copy.
// This removes the shared-mutable static config that BBRTestPGA used to hold.
//
// All tunables are reached through static Get*/Set* that route via Instance(),
// so callers automatically read/write their own thread's copy.

#include "globals.hh"
#include <iosfwd>

class BBRConfigMessenger;

class BBRConfigManager {
 public:
  static BBRConfigManager* Instance();   // thread-specific instance
  ~BBRConfigManager();                   // public for end-of-job cleanup

  // --- Emitter ---
  static G4double GetThermalT_K()  { return Instance()->fThermalT_K; }
  static void SetThermalT_K(G4double v);

  // --- Gun (values stored bare; PGA applies mm / eV) ---
  static G4bool   GetGunMode()     { return Instance()->fGunMode; }
  static G4double GetGunPosX_mm()  { return Instance()->fGunPosX_mm; }
  static G4double GetGunPosY_mm()  { return Instance()->fGunPosY_mm; }
  static G4double GetGunPosZ_mm()  { return Instance()->fGunPosZ_mm; }
  static G4double GetGunDirX()     { return Instance()->fGunDirX; }
  static G4double GetGunDirY()     { return Instance()->fGunDirY; }
  static G4double GetGunDirZ()     { return Instance()->fGunDirZ; }
  static G4double GetGunEnergy_eV(){ return Instance()->fGunEnergy_eV; }

  static void SetGunMode(G4bool v)      { Instance()->fGunMode = v; }
  static void SetGunPosX_mm(G4double v) { Instance()->fGunPosX_mm = v; }
  static void SetGunPosY_mm(G4double v) { Instance()->fGunPosY_mm = v; }
  static void SetGunPosZ_mm(G4double v) { Instance()->fGunPosZ_mm = v; }
  static void SetGunDirX(G4double v)    { Instance()->fGunDirX = v; }
  static void SetGunDirY(G4double v)    { Instance()->fGunDirY = v; }
  static void SetGunDirZ(G4double v)    { Instance()->fGunDirZ = v; }
  static void SetGunEnergy_eV(G4double v){ Instance()->fGunEnergy_eV = v; }

  // --- Detector (copper) ---
  static G4int    GetCuRRR()       { return Instance()->fCuRRR; }
  static G4double GetCuStageT_K()  { return Instance()->fCuStageT_K; }
  static void SetCuRRR(G4int rrr);
  static void SetCuStageT_K(G4double T_K);
  static void SetCuMaterial(const G4String& alias);  // named alias -> RRR + 4 K

  // --- Provenance ---
  static void Print(std::ostream& os) { Instance()->printConfig(os); }
  void printConfig(std::ostream& os) const;

 private:
  BBRConfigManager();                        // master/sequential: compiled defaults
  BBRConfigManager(const BBRConfigManager&); // clone from master (workers)

  BBRConfigManager(BBRConfigManager&&) = delete;
  BBRConfigManager& operator=(const BBRConfigManager&) = delete;
  BBRConfigManager& operator=(BBRConfigManager&&) = delete;

  G4double fThermalT_K;
  G4bool   fGunMode;
  G4double fGunPosX_mm, fGunPosY_mm, fGunPosZ_mm;
  G4double fGunDirX, fGunDirY, fGunDirZ;
  G4double fGunEnergy_eV;
  G4int    fCuRRR;
  G4double fCuStageT_K;

  BBRConfigMessenger* fMessenger;  // owned; one per thread instance
};

#endif
