#ifndef BBRPrimarySource_hh
#define BBRPrimarySource_hh

#include "ThermalSurface.hh"
#include "G4ParticleGun.hh"
#include <memory>

class G4Event;

// Shared source logic. Example-specific Geant4 user actions delegate to this
// class, while configuration remains per worker through BBRConfigManager.
class BBRPrimarySource {
 public:
  BBRPrimarySource();
  void Generate(G4Event* event);

 private:
  void BuildEmitter();

  std::unique_ptr<G4ParticleGun> fGun;
  ThermalSurface fSurface;
  G4ThreeVector fEmitterCenter_mm;
  G4ThreeVector fEmitterSize_mm;
};

#endif
