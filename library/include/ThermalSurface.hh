#ifndef ThermalSurface_h
#define ThermalSurface_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "BBEvt.hh"
#include "GetBBSpecCDF.hh"
#include "GeometricSurface.hh"
#include <vector>

class ThermalSurface {
 public:
  ThermalSurface();
  virtual ~ThermalSurface();

  // Add a box-shaped emitting surface.
  // center: world position; Wx/Wy/Wz: full extents (not half); in_out: true=outward;
  // rot1/rot2/rot3: rotations about the fixed world axes, applied in the order Z (rot1), then Y (rot2), then X (rot3); emissivity: 0–1.
  void AddBoxSurface(G4ThreeVector center,
                     G4double Wx, G4double Wy, G4double Wz,
                     G4bool in_out,
                     G4double rot1 = 0., G4double rot2 = 0., G4double rot3 = 0.,
                     G4double emissivity = 1.);

  void InitializeSpectrum(G4double temperature_K, G4double emin_eV, G4double emax_eV);
  G4double GetTemperature_K() const { return temp; }
  const GetBBSpecCDF& GetSpectrum() const { return BBSpecCDF; }
  G4double GetArea() const;
  G4double GetEffArea() const;

  // Remove all emitting surfaces (area totals reset). The Planck CDF and
  // temperature are untouched, so a caller can re-add surfaces and keep
  // sampling without re-initializing the spectrum.
  void ClearSurfaces();

  // Returns one BBEvt with energy (raw eV), position, and direction sampled
  // from the box surface. InitializeSpectrum must be called first.
  BBEvt GenEvt();

 private:
  GetBBSpecCDF BBSpecCDF;
  G4double temp = 0.;
  G4double area = 0.;
  G4double effArea = 0.;
  std::vector<GeometricSurface> surfaces;
};

#endif
