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
  std::unique_ptr<G4ParticleGun>     fGun;
  ThermalSurface                     fSurface;
  // static: shared across worker threads; messenger fires on master before beamOn
  static G4double                    fTemperature_K;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
#endif
