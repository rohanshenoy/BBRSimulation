#include "GetBBSpecCDF.hh"
#include "G4Exception.hh"
#include <cmath>

GetBBSpecCDF::GetBBSpecCDF()  {}
GetBBSpecCDF::~GetBBSpecCDF() {}

void GetBBSpecCDF::initialize(G4double temp, G4double emin, G4double emax)
{
  // Negated comparisons, so that a NaN argument is rejected too.
  if (!std::isfinite(temp) || !std::isfinite(emin) || !std::isfinite(emax) ||
      !(temp > 0.) || !(emin > 0.) || !(emin < emax))
    G4Exception("GetBBSpecCDF::initialize", "BBR017", FatalException,
                "needs finite T > 0 and 0 < emin < emax");
  // Build a complete new spectrum before replacing the previous valid one.
  std::vector<G4double> newX, newPDF, newCDF;
  newX.reserve(100001);
  newPDF.reserve(100001);
  newCDF.reserve(100001);

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
    newX.push_back(emin + G4double(count) * steps);
    if (!std::isfinite(newX[count]) ||
        (count > 0 && !(newX[count] > newX[count - 1])))
      G4Exception("GetBBSpecCDF::initialize", "BBR017", FatalException,
                  "Energy band is too narrow to make distinct CDF bins");
    // photon number spectrum: ∝ E²/(hc)²/(exp(E/kT)−1)
    G4double Bbody_y = 2. * newX[count] * newX[count] / hc / hc
                       / (std::exp(newX[count] / (k * temp)) - 1.);
    if (!std::isfinite(Bbody_y) || Bbody_y < 0.)
      G4Exception("GetBBSpecCDF::initialize", "BBR017", FatalException,
                  "Planck spectrum cannot be represented in the requested energy band");
    newPDF.push_back(Bbody_y);
    if (count > 0) sum += (newPDF[count] + newPDF[count - 1]) * steps / 2.;
    newCDF.push_back(sum);
    ++count;
  }

  if (!(sum > 0.) || !std::isfinite(sum))
    G4Exception("GetBBSpecCDF::initialize", "BBR017", FatalException,
                "Planck spectrum in the requested energy band cannot be normalized at this temperature");

  count = 0;
  while (count <= G4int(n_bin)) {
    newPDF[count] /= sum;
    newCDF[count] /= sum;
    if (!std::isfinite(newPDF[count]) || newPDF[count] < 0. ||
        !std::isfinite(newCDF[count]) || newCDF[count] < 0. ||
        newCDF[count] > 1. ||
        (count > 0 && newCDF[count] < newCDF[count - 1]))
      G4Exception("GetBBSpecCDF::initialize", "BBR017", FatalException,
                  "Normalized Planck spectrum is invalid");
    ++count;
  }
  x.swap(newX);
  pdf.swap(newPDF);
  cdf.swap(newCDF);
}
