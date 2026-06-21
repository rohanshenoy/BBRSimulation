#include "BBRTestPGA.hh"
#include "BBRConfigManager.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <cmath>

BBRTestPGA::BBRTestPGA()
  : G4VUserPrimaryGeneratorAction()
  , fGun(std::make_unique<G4ParticleGun>(1))
{
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhoton());

  // Planck emitter patch at x=-50mm.
  // Wx=1mm, Wy=Wz=20mm: ~91% of photons from x-faces, ~45% aimed at +x.
  fSurface.temp = BBRConfigManager::GetThermalT_K();
  fSurface.AddBoxSurface(
    G4ThreeVector(-50.*mm, 0., 0.),  // center
    1.*mm, 20.*mm, 20.*mm,           // Wx, Wy, Wz
    true,                            // in_out=1 -> outward emission
    0., 0., 0.,                      // no rotation
    1.0);                            // emissivity=1

  // Energy range 10 GHz-20 THz in eV (bare eV). 8.27e-2 eV = 20 THz.
  fSurface.BBSpecCDF.initialize(fSurface.temp, 4.14e-5, 8.27e-2);
}

void BBRTestPGA::GeneratePrimaries(G4Event* event)
{
  if (BBRConfigManager::GetGunMode()) {
    G4ThreeVector pos(BBRConfigManager::GetGunPosX_mm()*mm,
                      BBRConfigManager::GetGunPosY_mm()*mm,
                      BBRConfigManager::GetGunPosZ_mm()*mm);
    G4ThreeVector dir(BBRConfigManager::GetGunDirX(),
                      BBRConfigManager::GetGunDirY(),
                      BBRConfigManager::GetGunDirZ());
    dir = dir.unit();
    G4ThreeVector perp = dir.orthogonal().unit();
    G4double phi = G4UniformRand() * CLHEP::twopi;
    G4ThreeVector pol = std::cos(phi)*perp + std::sin(phi)*dir.cross(perp).unit();
    fGun->SetParticlePosition(pos);
    fGun->SetParticleMomentumDirection(dir);
    fGun->SetParticleEnergy(BBRConfigManager::GetGunEnergy_eV() * eV);
    fGun->SetParticlePolarization(pol);
    fGun->GeneratePrimaryVertex(event);
    return;
  }

  // Re-initialize CDF if the emitter temperature changed via messenger.
  const G4double T = BBRConfigManager::GetThermalT_K();
  if (fSurface.temp != T) {
    fSurface.temp = T;
    fSurface.BBSpecCDF.initialize(T, 4.14e-5, 8.27e-2);
  }

  BBEvt evt = fSurface.GenEvt();

  G4ThreeVector dir  = evt.direction.unit();
  G4ThreeVector perp = dir.orthogonal().unit();
  G4double      phi  = G4UniformRand() * CLHEP::twopi;
  G4ThreeVector pol  = std::cos(phi) * perp
                     + std::sin(phi) * dir.cross(perp).unit();

  fGun->SetParticlePosition(evt.position);
  fGun->SetParticleMomentumDirection(dir);
  fGun->SetParticleEnergy(evt.energy * eV);  // raw eV -> Geant4 internal units
  fGun->SetParticlePolarization(pol);
  fGun->GeneratePrimaryVertex(event);
}
