#include "GetBBSpecCDF.hh"
#include <cmath>

GetBBSpecCDF::GetBBSpecCDF()  {}
GetBBSpecCDF::~GetBBSpecCDF() {}

void GetBBSpecCDF::initialize(G4double temp, G4double emin, G4double emax)
{
  // Clear existing data before re-initializing (supports runtime temperature changes).
  x.clear();
  pdf.clear();
  cdf.clear();

  const G4double n_bin  = 100000;
  const G4double erange = emax - emin;
  const G4double steps  = erange / n_bin;

  // Constants in eV/SI — consistent with x[] in raw eV.
  const G4double k  = 8.6173e-5;   // Boltzmann, eV/K
  const G4double h  = 4.1357e-15;  // Planck, eV·s
  const G4double c  = 2.9979e8;    // m/s
  const G4double hc = h * c;

  G4int    count = 0;
  G4double sum   = 0.;

  while (count <= G4int(n_bin)) {
    x.push_back(emin + G4double(count) * steps);
    // photon number spectrum: ∝ E²/(hc)²/(exp(E/kT)−1)
    G4double Bbody_y = 2. * x[count] * x[count] / hc / hc
                       / (std::exp(x[count] / (k * temp)) - 1.);
    pdf.push_back(Bbody_y);
    if (count > 0) sum += (pdf[count] + pdf[count - 1]) * steps / 2.;
    cdf.push_back(sum);
    ++count;
  }

  count = 0;
  while (count <= G4int(n_bin)) {
    pdf[count] /= sum;
    cdf[count] /= sum;
    ++count;
  }
}
