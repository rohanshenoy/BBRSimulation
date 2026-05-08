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
  // rot1/rot2/rot3: Euler angles (rotZ, rotY', rotX''); emissivity: 0–1.
  void AddBoxSurface(G4ThreeVector center,
                     G4double Wx, G4double Wy, G4double Wz,
                     G4bool in_out,
                     G4double rot1 = 0., G4double rot2 = 0., G4double rot3 = 0.,
                     G4double emissivity = 1.);

  GetBBSpecCDF BBSpecCDF;   // call BBSpecCDF.initialize(T, emin, emax) before use
  G4double     temp = 0.;  // K — informational
  G4double     area    = 0.;
  G4double     effArea = 0.;
  G4double     GetArea();
  G4double     GetEffArea();

  // Returns one BBEvt with energy (raw eV), position, and direction sampled
  // from the box surface. BBSpecCDF must be initialized before calling.
  BBEvt GenEvt();

 private:
  std::vector<GeometricSurface> surfaces;
};

#endif
