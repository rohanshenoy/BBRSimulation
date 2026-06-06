#ifndef BBRTestDetectorConstruction_hh
#define BBRTestDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include <memory>

// Test geometry: 50 cm world, 4 mm Cu slab at x=2 mm with two vacuum_wg crack daughters.
// The Cu wall material is built from (fCuRRR, fStageT_K) at Construct() time via
// BBRMaterials::GetCopper(RRR, T_K).  Both can be set before /run/initialize.
//
// Messenger commands — all must be issued before /run/initialize:
//
//   /bbr/det/setCuMaterial <alias>
//     Named shortcut.  Resolves to an RRR value; does not accept raw numbers.
//     Valid aliases:  OFHC_Cu (RRR=100)  |  OF_Cu (RRR=3)  |  HP_Cu (RRR=6)
//
//   /bbr/det/setCuRRR <N>
//     Direct integer RRR (Residual Resistance Ratio, >= 1).
//     Physics: sigma_DC = RRR * sigma_RT, where sigma_RT = 5.96e7 S/m for all Cu
//     at 273 K (universal constant — users never need to supply it).
//     Examples: 50 (annealed OFHC), 100 (standard OFHC), 300 (ultra-pure crystal).
//
//   /bbr/det/setCuStageT <T> K
//     Temperature of the cryostat stage this Cu component belongs to.
//     Default: 4 K.  Affects reflectance table via Drude scattering time tau.
//     Below 50 K: sigma_DC = RRR * sigma_RT (impurity term dominates; phonons frozen).
//     Above 50 K: Matthiessen's rule with sigma_phonon ~ sigma_RT * 273/T.
//     Note: above ~50 K the 1/T approximation breaks down (umklapp/Bloch-Grüneisen);
//     for shield layers at 10–50 K the Matthiessen term is approximate.
class BBRTestDetectorConstruction : public G4VUserDetectorConstruction {
public:
  BBRTestDetectorConstruction();
  ~BBRTestDetectorConstruction() override = default;
  G4VPhysicalVolume* Construct() override;

  void SetCuMaterial(const G4String& alias);  // named alias → sets fCuRRR
  void SetCuRRR(G4int rrr);                   // direct RRR; validates >= 1
  void SetCuStageT(G4double T_K);             // temperature stage [K]

private:
  G4int    fCuRRR    = 100;  // Residual Resistance Ratio; default OFHC_Cu
  G4double fStageT_K = 4.0;  // temperature [K]; default cryogenic baseline
  std::unique_ptr<G4GenericMessenger> fMessenger;
};

#endif
