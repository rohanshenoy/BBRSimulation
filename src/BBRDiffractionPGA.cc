#include "BBRDiffractionPGA.hh"

#include "G4GenericMessenger.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"

BBRDiffractionPGA::BBRDiffractionPGA(G4double gunZ_mm)
{
  fGun = new G4ParticleGun(1);
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
  fGun->SetParticlePosition(G4ThreeVector(-20.*mm, 0., gunZ_mm));
  fGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fGun->SetParticleEnergy(2.068e-3 * eV);  // 500 GHz

  // Polarization: phat = (0, -1/√2, -1/√2).
  // At normal incidence: E_theta = phat·ê_θ_in = 1/√2, E_phi = 1/√2.
  // Expected: T = 0.5*T_Ephi0 + 0.5*T_Ephi1 ≈ 0.5*1.055 + 0 ≈ 0.528.
  fGun->SetParticlePolarization(G4ThreeVector(0., -M_SQRT1_2, -M_SQRT1_2));

  fMessenger = new G4GenericMessenger(this, "/bbr/gun/", "BBR gun controls");
  fMessenger->DeclareMethodWithUnit("setZ", "mm", &BBRDiffractionPGA::SetGunZ,
                                    "Set gun Z position in mm (selects which crack to hit)");
}

BBRDiffractionPGA::~BBRDiffractionPGA()
{
  delete fMessenger;
  delete fGun;
}

void BBRDiffractionPGA::SetGunZ(G4double z)
{
  fGun->SetParticlePosition(G4ThreeVector(-20.*mm, 0., z));
}

void BBRDiffractionPGA::GeneratePrimaries(G4Event* event)
{
  fGun->GeneratePrimaryVertex(event);
}
