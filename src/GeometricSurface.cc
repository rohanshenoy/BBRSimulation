#include "GeometricSurface.hh"
#include <cmath>

GeometricSurface::GeometricSurface()  {}
GeometricSurface::~GeometricSurface() {}

G4double GeometricSurface::CalculateArea()
{
  const G4double PI = 3.14159265358979323846;
  G4double A = 0.;
  switch (type) {
    case 1:  // tube
      A = 2. * PI * radius * (height * 2.);
      if (lid) A += PI * radius * radius * 2.;
      break;
    case 2:  // disc
      A = PI * (radius * radius - r_in * r_in);
      if (both_side) A *= 2.;
      break;
    case 3:  // box
      A = 2. * (Wx * Wy + Wy * Wz + Wz * Wx);
      break;
    case 4:  // sphere
      A = radius * radius * (std::cos(theta1) - std::cos(theta2)) * (phi2 - phi1);
      break;
  }
  area = A;
  return A;
}
