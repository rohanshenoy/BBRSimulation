#ifndef TestWorldPrimaryGeneratorAction_hh
#define TestWorldPrimaryGeneratorAction_hh

#include "ThermalSurface.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <memory>

// Dual-mode primary generator. Reads its settings (emitter T and box
// geometry, gun mode/pos/dir/energy) from BBRConfigManager each event; owns no
// config state. The Planck emitter box is rebuilt whenever
// /bbr/thermal/emitterCenter or /bbr/thermal/emitterSize changes.
class TestWorldPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  TestWorldPrimaryGeneratorAction();
  ~TestWorldPrimaryGeneratorAction() override = default;
  void GeneratePrimaries(G4Event*) override;

 private:
  void BuildEmitter();   // (re)build fSurface's box from BBRConfigManager

  std::unique_ptr<G4ParticleGun> fGun;
  ThermalSurface                 fSurface;
  G4ThreeVector                  fEmitterCenter_mm;   // geometry fSurface was built with
  G4ThreeVector                  fEmitterSize_mm;
};
#endif
