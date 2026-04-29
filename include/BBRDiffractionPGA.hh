#ifndef BBRDiffractionPGA_hh
#define BBRDiffractionPGA_hh

#include "G4Types.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// Primary generator for the diffraction smoke test.
// Fires one opticalphoton per event:
//   position  : (-20, 0, z) mm  — upstream of the TEM_waveguide_crack entry face
//   direction : (+1,  0, 0)      — normal incidence
//   energy    : 2.068 meV        — 500 GHz
//   polaris.  : (0, -1/√2, -1/√2) — equal E_theta / E_phi at normal incidence (T ≈ 50%)
// Gun Z can be changed at runtime via /bbr/gun/setZ <value> mm.
class G4Event;
class G4GenericMessenger;
class G4ParticleGun;

class BBRDiffractionPGA : public G4VUserPrimaryGeneratorAction
{
 public:
  explicit BBRDiffractionPGA(G4double gunZ_mm = 0.);
  ~BBRDiffractionPGA() override;

  void GeneratePrimaries(G4Event* event) override;
  void SetGunZ(G4double z);

 private:
  G4ParticleGun*      fGun       = nullptr;
  G4GenericMessenger* fMessenger = nullptr;
};

#endif
