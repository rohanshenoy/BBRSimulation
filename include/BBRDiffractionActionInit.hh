#ifndef BBRDiffractionActionInit_hh
#define BBRDiffractionActionInit_hh

#include "G4Types.hh"
#include "G4VUserActionInitialization.hh"

// Minimal action initialization for the diffraction smoke test.
// Uses BBRDiffractionPGA (normal-incidence opticalphoton) and
// RunAction(nullptr); skips OpNovice2 SteppingAction to avoid the
// group-velocity assertion that fires when both volumes share the same
// G4_Galactic material.
// gunZ_mm: world-frame z position of the gun (0 = crack1, 3 mm = crack2).
class BBRDiffractionActionInit : public G4VUserActionInitialization
{
 public:
  explicit BBRDiffractionActionInit(G4double gunZ_mm = 0.);
  ~BBRDiffractionActionInit() override = default;

  void BuildForMaster() const override;
  void Build()          const override;

 private:
  G4double fGunZ_mm = 0.;
};

#endif
