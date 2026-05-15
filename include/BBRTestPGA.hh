#ifndef BBRTestPGA_hh
#define BBRTestPGA_hh

#include "ThermalSurface.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GenericMessenger.hh"
#include <memory>

class BBRTestPGA : public G4VUserPrimaryGeneratorAction {
public:
  BBRTestPGA();
  ~BBRTestPGA() override = default;
  void GeneratePrimaries(G4Event*) override;

private:
  std::unique_ptr<G4ParticleGun>      fGun;
  ThermalSurface                      fSurface;

  // static: shared across worker threads; messengers fire on master before beamOn
  static G4double fTemperature_K;

  // Gun mode: fire a fixed-position/direction photon instead of Planck sampling.
  // Useful for shooting directly at a crack to get angular statistics.
  static G4bool   fGunMode;
  static G4double fGunPosX_mm, fGunPosY_mm, fGunPosZ_mm;
  static G4double fGunDirX, fGunDirY, fGunDirZ;
  static G4double fGunEnergy_eV;   // default 2.07e-3 eV = 500 GHz

  std::unique_ptr<G4GenericMessenger> fMessenger;
  std::unique_ptr<G4GenericMessenger> fGunMessenger;
};
#endif
