#ifndef GetBBSpecCDF_h
#define GetBBSpecCDF_h 1

#include "globals.hh"
#include <vector>

// Builds the CDF for the Planck photon-number spectrum B ∝ ν²/(exp(hν/kT)-1).
// All energies are in raw eV (NOT Geant4 internal units).
// Call initialize(T_K, emin_eV, emax_eV) once before first use.
class GetBBSpecCDF {
 public:
  GetBBSpecCDF();
  virtual ~GetBBSpecCDF();

  void initialize(G4double temperature_K, G4double emin_eV, G4double emax_eV);

  std::vector<G4double> x;    // energy axis, raw eV
  std::vector<G4double> pdf;  // normalised PDF
  std::vector<G4double> cdf;  // normalised CDF
};

#endif
