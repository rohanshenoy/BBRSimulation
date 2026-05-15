#include "BBRTestPGA.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <cmath>

// Static member: shared across all worker-thread PGA instances so that
// /bbr/thermal/setT (which updates the master-thread instance) propagates
// to all workers before /run/beamOn starts.
G4double BBRTestPGA::fTemperature_K = 4.0;

BBRTestPGA::BBRTestPGA()
  : G4VUserPrimaryGeneratorAction()
  , fGun(std::make_unique<G4ParticleGun>(1))
{
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhoton());

  // Planck 4K emitter patch at x=-50mm.
  // Wx=1mm, Wy=Wz=20mm: ~91% of photons from x-faces, ~45% aimed at +x toward Cu wall.
  fSurface.temp = 4.0;
  fSurface.AddBoxSurface(
    G4ThreeVector(-50.*mm, 0., 0.),  // center
    1.*mm, 20.*mm, 20.*mm,           // Wx, Wy, Wz
    true,                            // in_out=1 → outward emission
    0., 0., 0.,                      // no rotation
    1.0);                            // emissivity=1

  // Energy range 50 GHz–20 THz in eV (bare eV, not Geant4 internal units)
  fSurface.BBSpecCDF.initialize(fSurface.temp, 4.14e-5, 2.07e-2);

  fMessenger = std::make_unique<G4GenericMessenger>(this, "/bbr/thermal/", "Thermal emitter settings");
  fMessenger->DeclareProperty("setT", fTemperature_K, "Emitter temperature [K]")
            .SetParameterName("T", false)
            .SetRange("T > 0")
            .SetDefaultValue("4.0");
}

void BBRTestPGA::GeneratePrimaries(G4Event* event)
{
  // Re-initialize CDF if temperature changed via messenger
  if (fSurface.temp != fTemperature_K) {
    fSurface.temp = fTemperature_K;
    fSurface.BBSpecCDF.initialize(fTemperature_K, 4.14e-5, 2.07e-2);
  }

  BBEvt evt = fSurface.GenEvt();

  G4ThreeVector dir  = evt.direction.unit();
  G4ThreeVector perp = dir.orthogonal().unit();
  G4double      phi  = G4UniformRand() * CLHEP::twopi;
  G4ThreeVector pol  = std::cos(phi) * perp
                     + std::sin(phi) * dir.cross(perp).unit();

  fGun->SetParticlePosition(evt.position);
  fGun->SetParticleMomentumDirection(dir);
  fGun->SetParticleEnergy(evt.energy * eV);  // raw eV → Geant4 internal units
  fGun->SetParticlePolarization(pol);
  fGun->GeneratePrimaryVertex(event);
}
