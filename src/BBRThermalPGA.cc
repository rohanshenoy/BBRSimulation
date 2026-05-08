#include "BBRThermalPGA.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>

BBRThermalPGA::BBRThermalPGA()
  : fGun(new G4ParticleGun(1))
{
  fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());

  // 4K patch: 1cm × 1cm × 1mm box at origin, emitting outward (in_out=true).
  fSurface.temp = 4.;
  fSurface.AddBoxSurface(
    G4ThreeVector(0., 0., 0.),   // center (world coords, Geant4 units)
    1.*cm, 1.*cm, 1.*mm,          // Wx, Wy, Wz (Geant4 units)
    true                          // in_out=true → outward emission
  );
  // Energy range: 10 GHz – 5 THz. Passed as bare eV numbers (not * eV) because
  // GetBBSpecCDF::initialize uses eV-based constants internally.
  const G4double emin_eV = 4.14e-5;  // 10 GHz in eV
  const G4double emax_eV = 2.07e-2;  // 5 THz in eV
  fSurface.BBSpecCDF.initialize(fSurface.temp, emin_eV, emax_eV);

  fOut.open("planck_output.csv");
  fOut << "energy_eV\n";
}

BBRThermalPGA::~BBRThermalPGA() { delete fGun; }

void BBRThermalPGA::GeneratePrimaries(G4Event* event)
{
  BBEvt evt = fSurface.GenEvt();

  // Random linear polarization perpendicular to direction.
  // Required for G4OpticalPhoton — YYC's code does not set this.
  G4ThreeVector perpA = evt.direction.orthogonal().unit();
  G4ThreeVector perpB = evt.direction.cross(perpA).unit();
  G4double psi = CLHEP::twopi * G4UniformRand();
  G4ThreeVector pol = std::cos(psi) * perpA + std::sin(psi) * perpB;

  // evt.energy is in raw eV; convert to Geant4 internal units (* eV).
  fGun->SetParticleEnergy(evt.energy * eV);
  fGun->SetParticlePosition(evt.position);
  fGun->SetParticleMomentumDirection(evt.direction);
  fGun->SetParticlePolarization(pol);
  fGun->GeneratePrimaryVertex(event);

  // CSV stores raw eV (what we read back in Python).
  fOut << evt.energy << "\n";
}
