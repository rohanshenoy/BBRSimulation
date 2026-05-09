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

  G4int         type;      // 1=tube, 2=disc, 3=box, 4=sphere
  G4ThreeVector position;
  G4ThreeVector direction;
  G4double      area;
  G4double      CalculateArea();

  // tube fields
  G4ThreeVector center;
  G4double      radius;
  G4double      height;
  G4bool        in_out;
  G4bool        lid;
  G4double      rot1;
  G4double      rot2;
  G4double      rot3;

  // disc fields
  G4double      r_in;
  G4bool        both_side;

  // box fields
  G4double      Wx;
  G4double      Wy;
  G4double      Wz;

  // sphere fields
  G4double      theta1;
  G4double      theta2;
  G4double      phi1;
  G4double      phi2;
};

#endif
