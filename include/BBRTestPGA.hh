#ifndef BBRTestPGA_hh
#define BBRTestPGA_hh

#include "ThermalSurface.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <memory>

// Dual-mode primary generator. Reads its settings (emitter T, gun mode/pos/
// dir/energy) from BBRConfigManager each event; owns no config state.
class BBRTestPGA : public G4VUserPrimaryGeneratorAction {
 public:
  BBRTestPGA();
  ~BBRTestPGA() override = default;
  void GeneratePrimaries(G4Event*) override;

 private:
  std::unique_ptr<G4ParticleGun> fGun;
  ThermalSurface                 fSurface;
};
#endif
