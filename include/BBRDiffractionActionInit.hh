#ifndef BBRDiffractionActionInit_hh
#define BBRDiffractionActionInit_hh

#include "G4VUserActionInitialization.hh"

// Minimal action initialization for the diffraction smoke test.
// Uses BBRDiffractionPGA (hardcoded normal-incidence opticalphoton) and
// RunAction(nullptr); skips OpNovice2 SteppingAction to avoid the
// group-velocity assertion that fires when both volumes share the same
// G4_Galactic material.
class BBRDiffractionActionInit : public G4VUserActionInitialization
{
 public:
  BBRDiffractionActionInit()  = default;
  ~BBRDiffractionActionInit() override = default;

  void BuildForMaster() const override;
  void Build()          const override;
};

#endif
