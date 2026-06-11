#ifndef GeometricSurface_h
#define GeometricSurface_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

// Holds geometry parameters for one primitive surface element.
// YYC's GeometricSurface with CDMSSnolabIRBackgroundDetectorConstruction
// inheritance removed (that class is not present and its methods are unused).
class GeometricSurface {
 public:
  GeometricSurface();
  virtual ~GeometricSurface();

  G4int         type = 0;      // 1=tube, 2=disc, 3=box, 4=sphere
  G4ThreeVector position;
  G4ThreeVector direction;
  G4double      area = 0.;
  G4double      emissivity = 1.;   // 0–1, weights surface selection in GenEvt
  G4double      CalculateArea();

  // tube fields
  G4ThreeVector center;
  G4double      radius = 0.;
  G4double      height = 0.;
  G4bool        in_out = false;
  G4bool        lid = false;
  G4double      rot1 = 0.;
  G4double      rot2 = 0.;
  G4double      rot3 = 0.;

  // disc fields
  G4double      r_in = 0.;
  G4bool        both_side = false;

  // box fields
  G4double      Wx = 0.;
  G4double      Wy = 0.;
  G4double      Wz = 0.;

  // sphere fields
  G4double      theta1 = 0.;
  G4double      theta2 = 0.;
  G4double      phi1 = 0.;
  G4double      phi2 = 0.;
};

#endif
