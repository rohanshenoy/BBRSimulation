#ifndef BBRThermalPGA_hh
#define BBRThermalPGA_hh

#include "ThermalSurface.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include <fstream>

class G4Event;

class BBRThermalPGA : public G4VUserPrimaryGeneratorAction {
 public:
  BBRThermalPGA();
  ~BBRThermalPGA() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun;
  ThermalSurface fSurface;
  std::ofstream  fOut;    // writes planck_output.csv
};

#endif
