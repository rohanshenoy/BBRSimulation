#include "TestWorldPrimaryGeneratorAction.hh"
#include "BBRConfigManager.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <atomic>
#include <cmath>

namespace {
// The Planck CDF band is fixed at 10 GHz-20 THz (4.14e-5 to 8.27e-2 eV).
// Outside 2.2-117 K more than 1 % of the photon-number spectrum lies outside
// it and is never emitted (13.8 % at 0.5 K, 3.9 % at 150 K). Once per process.
void WarnIfOutOfBand(G4double T_K, const char* origin)
{
  static std::atomic<bool> warned{false};
  if ((T_K < 2.2 || T_K > 117.) && !warned.exchange(true)) {
    G4ExceptionDescription ed;
    ed << "Emitter temperature " << T_K << " K is outside 2.2-117 K: more than 1 % of the "
       << "Planck photon-number spectrum lies outside the 10 GHz-20 THz band and is not emitted.";
    G4Exception(origin, "BBR021", JustWarning, ed);
  }
}
}  // namespace

TestWorldPrimaryGeneratorAction::TestWorldPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
  , fGun(std::make_unique<G4ParticleGun>(1))
{
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhoton());

  // Planck emitter box from BBRConfigManager. Test-world default: a 1x20x20 mm
  // patch centred at x=-50 mm (~91% of photons from the x-faces, ~45% aimed
  // at +x). Other geometries (e.g. the light pipe) set their own box via
  // /bbr/thermal/emitterCenter and /bbr/thermal/emitterSize.
  fSurface.temp = BBRConfigManager::GetThermalT_K();
  BuildEmitter();

  // Energy range 10 GHz-20 THz in eV (bare eV). 8.27e-2 eV = 20 THz.
  fSurface.BBSpecCDF.initialize(fSurface.temp, 4.14e-5, 8.27e-2);
}

void TestWorldPrimaryGeneratorAction::BuildEmitter()
{
  fEmitterCenter_mm = BBRConfigManager::GetEmitterCenter_mm();
  fEmitterSize_mm   = BBRConfigManager::GetEmitterSize_mm();
  fSurface.ClearSurfaces();
  fSurface.AddBoxSurface(
    fEmitterCenter_mm * mm,                         // center
    fEmitterSize_mm.x() * mm,                       // Wx
    fEmitterSize_mm.y() * mm,                       // Wy
    fEmitterSize_mm.z() * mm,                       // Wz
    true,                                           // in_out=1 -> outward emission
    0., 0., 0.,                                     // no rotation
    1.0);                                           // emissivity=1
}

void TestWorldPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (BBRConfigManager::GetGunMode()) {
    G4ThreeVector pos(BBRConfigManager::GetGunPosX_mm()*mm,
                      BBRConfigManager::GetGunPosY_mm()*mm,
                      BBRConfigManager::GetGunPosZ_mm()*mm);
    G4ThreeVector dir(BBRConfigManager::GetGunDirX(),
                      BBRConfigManager::GetGunDirY(),
                      BBRConfigManager::GetGunDirZ());
    dir = dir.unit();
    // Polarization: a fixed vector from /bbr/gun/pol (projected perpendicular
    // to the direction and normalised), or, for the zero vector (default),
    // random in the plane perpendicular to the direction.
    G4ThreeVector pol;
    const G4ThreeVector polReq  = BBRConfigManager::GetGunPol();
    const G4ThreeVector polPerp = polReq - polReq.dot(dir) * dir;
    if (polReq.mag2() > 0. && polPerp.mag() > 1e-9) {
      pol = polPerp.unit();
    } else {
      if (polReq.mag2() > 0.) {
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true)) {
          G4Exception("TestWorldPrimaryGeneratorAction::GeneratePrimaries", "BBR010", JustWarning,
                      "/bbr/gun/pol is parallel to the gun direction; "
                      "using random polarization.");
        }
      }
      G4ThreeVector perp = dir.orthogonal().unit();
      G4double phi = G4UniformRand() * CLHEP::twopi;
      pol = std::cos(phi)*perp + std::sin(phi)*dir.cross(perp).unit();
    }
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
  WarnIfOutOfBand(T, "TestWorldPrimaryGeneratorAction::GeneratePrimaries");
  // Rebuild the emitter box if its geometry changed via messenger.
  if (BBRConfigManager::GetEmitterCenter_mm() != fEmitterCenter_mm ||
      BBRConfigManager::GetEmitterSize_mm()   != fEmitterSize_mm) {
    BuildEmitter();
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
