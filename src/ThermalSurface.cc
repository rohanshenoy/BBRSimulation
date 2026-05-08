#include "ThermalSurface.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

ThermalSurface::ThermalSurface()  {}
ThermalSurface::~ThermalSurface() {}

void ThermalSurface::AddBoxSurface(G4ThreeVector center,
                                    G4double Wx, G4double Wy, G4double Wz,
                                    G4bool in_out,
                                    G4double rot1, G4double rot2, G4double rot3,
                                    G4double emissivity)
{
  GeometricSurface s;
  s.type   = 3;
  s.center = center;
  s.Wx     = Wx;
  s.Wy     = Wy;
  s.Wz     = Wz;
  s.in_out = in_out;
  s.rot1   = rot1;
  s.rot2   = rot2;
  s.rot3   = rot3;
  s.CalculateArea();

  surfaces.push_back(s);
  area    += s.area;
  effArea += s.area * emissivity;
}

G4double ThermalSurface::GetArea()    { return area; }
G4double ThermalSurface::GetEffArea() { return effArea; }

BBEvt ThermalSurface::GenEvt()
{
  G4ThreeVector X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1);
  BBEvt thisEvt;

  // ---- energy: inverse CDF with quadratic interpolation (YYC exact) ----
  G4double BBSpecProbBelow = G4UniformRand();
  G4int    nCDF            = (G4int)BBSpecCDF.cdf.size();
  int      idx             = nCDF - 2;
  for (int i = 0; i < nCDF; ++i) {
    if (BBSpecCDF.cdf[i] > BBSpecProbBelow) { idx = i - 1; break; }
  }
  if (idx < 0)        idx = 0;
  if (idx >= nCDF - 1) idx = nCDF - 2;

  G4double a = (BBSpecCDF.pdf[idx + 1] - BBSpecCDF.pdf[idx])
               / (BBSpecCDF.x[idx + 1] - BBSpecCDF.x[idx]) / 2.;
  G4double b = BBSpecCDF.pdf[idx];
  G4double c = BBSpecCDF.cdf[idx] - BBSpecProbBelow;

  if (a != 0.)
    thisEvt.energy = BBSpecCDF.x[idx] + (-b + std::sqrt(b * b - 4. * a * c)) / 2. / a;
  else
    thisEvt.energy = BBSpecCDF.x[idx] - c / b;

  // ---- surface selection: area-weighted ----
  if (surfaces.empty()) {
    G4Exception("ThermalSurface::GenEvt", "BBR001", FatalException,
                "No surfaces added. Call AddBoxSurface before GenEvt.");
  }
  G4double totalArea  = 0.;
  int      N_surfaces = (int)surfaces.size();
  for (int i = 0; i < N_surfaces; ++i) totalArea += surfaces[i].area;

  G4double prob = G4UniformRand() * totalArea;
  int idx_s = 0;
  for (int i = 0; i < N_surfaces; ++i) {
    prob -= surfaces[i].area;
    if (prob <= 0.) { idx_s = i; break; }
  }

  // ---- box emission: port of YYC case 3 ----
  G4double      Wx         = surfaces[idx_s].Wx;
  G4double      Wy         = surfaces[idx_s].Wy;
  G4double      Wz         = surfaces[idx_s].Wz;
  G4bool        Box_in_out = surfaces[idx_s].in_out;
  G4double      rot1       = surfaces[idx_s].rot1;
  G4double      rot2       = surfaces[idx_s].rot2;
  G4double      rot3       = surfaces[idx_s].rot3;
  G4ThreeVector BoxCenter  = surfaces[idx_s].center;

  G4double tot_area     = 2. * (Wx * Wy + Wy * Wz + Wz * Wx);
  G4double face_pn      = (G4UniformRand() > 0.5) ? 1. : -1.;
  G4double radiation_pn = -face_pn;                       // default: inward
  if (Box_in_out == 1) radiation_pn = face_pn;            // in_out=1 → outward

  G4double xx, yy, zz;
  G4double polTheta  = G4UniformRand() * 90. * deg;      // YYC: uniform in θ
  G4double polPhi    = G4UniformRand() * 360. * deg;
  G4double face_rand = G4UniformRand();

  if (face_rand < 2. * Wy * Wz / tot_area) {
    // x-face
    xx = Wx / 2. * face_pn;
    yy = (G4UniformRand() - 0.5) * Wy;
    zz = (G4UniformRand() - 0.5) * Wz;
    thisEvt.position  = G4ThreeVector(xx, yy, zz)
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X) + BoxCenter;
    thisEvt.direction = G4ThreeVector(std::cos(polTheta) * radiation_pn,
                                      std::sin(polTheta) * std::cos(polPhi),
                                      std::sin(polTheta) * std::sin(polPhi))
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X);
  } else if (face_rand < 2. * (Wy * Wz + Wz * Wx) / tot_area) {
    // y-face
    xx = (G4UniformRand() - 0.5) * Wx;
    yy = Wy / 2. * face_pn;
    zz = (G4UniformRand() - 0.5) * Wz;
    thisEvt.position  = G4ThreeVector(xx, yy, zz)
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X) + BoxCenter;
    thisEvt.direction = G4ThreeVector(std::sin(polTheta) * std::cos(polPhi),
                                      std::cos(polTheta) * radiation_pn,
                                      std::sin(polTheta) * std::sin(polPhi))
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X);
  } else {
    // z-face
    xx = (G4UniformRand() - 0.5) * Wx;
    yy = (G4UniformRand() - 0.5) * Wy;
    zz = Wz / 2. * face_pn;
    thisEvt.position  = G4ThreeVector(xx, yy, zz)
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X) + BoxCenter;
    thisEvt.direction = G4ThreeVector(std::sin(polTheta) * std::cos(polPhi),
                                      std::sin(polTheta) * std::sin(polPhi),
                                      std::cos(polTheta) * radiation_pn)
                          .rotate(rot1, Z).rotate(rot2, Y).rotate(rot3, X);
  }

  return thisEvt;
}
