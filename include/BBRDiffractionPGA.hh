#ifndef BBRDiffractionPGA_hh
#define BBRDiffractionPGA_hh

#include "G4Types.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// Primary generator for the diffraction smoke test.
// Fires one opticalphoton per event:
//   position  : (-20, 0, 0) mm  — upstream of the TEM_waveguide_crack entry face at x=-0.5 mm
//   direction : (+1,  0, 0)      — normal incidence
//   energy    : 2.068 meV        — 500 GHz
//   polaris.  : (0, -1/√2, -1/√2) — equal E_theta / E_phi at normal incidence (T ≈ 50%)
class G4Event;
class G4ParticleGun;

class BBRDiffractionPGA : public G4VUserPrimaryGeneratorAction
{
 public:
  explicit BBRDiffractionPGA(G4double gunZ_mm = 0.);
  ~BBRDiffractionPGA() override;

  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun = nullptr;
};

#endif
