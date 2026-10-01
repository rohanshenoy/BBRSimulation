#include "ThermalSurface.hh"
#include "globals.hh"
#include "G4Exception.hh"
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
  if (!std::isfinite(center.x()) || !std::isfinite(center.y()) ||
      !std::isfinite(center.z()) || !std::isfinite(Wx) ||
      !std::isfinite(Wy) || !std::isfinite(Wz) ||
      !(Wx > 0.) || !(Wy > 0.) || !(Wz > 0.) ||
      !std::isfinite(rot1) || !std::isfinite(rot2) || !std::isfinite(rot3) ||
      !std::isfinite(emissivity) || emissivity < 0. || emissivity > 1.)
    G4Exception("ThermalSurface::AddBoxSurface", "BBR023", FatalException,
                "Box needs finite center/rotation, positive extents and emissivity in [0,1]");
  GeometricSurface s;
  s.type       = 3;
  s.center     = center;
  s.Wx         = Wx;
  s.Wy         = Wy;
  s.Wz         = Wz;
  s.in_out     = in_out;
  s.rot1       = rot1;
  s.rot2       = rot2;
  s.rot3       = rot3;
  s.emissivity = emissivity;
  s.CalculateArea();
  if (!std::isfinite(s.area) || !std::isfinite(area + s.area) ||
      !std::isfinite(effArea + s.area * emissivity))
    G4Exception("ThermalSurface::AddBoxSurface", "BBR023", FatalException,
                "Box area or emitting area is too large");

  surfaces.push_back(s);
  area    += s.area;
  effArea += s.area * emissivity;
}

void ThermalSurface::InitializeSpectrum(G4double temperature_K, G4double emin_eV, G4double emax_eV)
{
  BBSpecCDF.initialize(temperature_K, emin_eV, emax_eV);
  temp = temperature_K;
}

G4double ThermalSurface::GetArea() const    { return area; }
G4double ThermalSurface::GetEffArea() const { return effArea; }

void ThermalSurface::ClearSurfaces()
{
  surfaces.clear();
  area    = 0.;
  effArea = 0.;
}

BBEvt ThermalSurface::GenEvt()
{
  if (BBSpecCDF.CDF().size() < 2) G4Exception("ThermalSurface::GenEvt", "BBR018", FatalException, "Planck CDF not initialized: call InitializeSpectrum first");
  G4ThreeVector X(1, 0, 0), Y(0, 1, 0), Z(0, 0, 1);
  BBEvt thisEvt;

  // ---- energy: inverse CDF with quadratic interpolation (YYC exact) ----
  G4double BBSpecProbBelow = G4UniformRand();
  const auto& cdf = BBSpecCDF.CDF();
  const auto& pdf = BBSpecCDF.PDF();
  const auto& x = BBSpecCDF.EnergyAxis();
  G4int    nCDF            = (G4int)cdf.size();
  int      idx             = nCDF - 2;
  for (int i = 0; i < nCDF; ++i) {
    if (cdf[i] > BBSpecProbBelow) { idx = i - 1; break; }
  }
  if (idx < 0)        idx = 0;
  if (idx >= nCDF - 1) idx = nCDF - 2;

  G4double a = (pdf[idx + 1] - pdf[idx])
               / (x[idx + 1] - x[idx]) / 2.;
  G4double b = pdf[idx];
  G4double c = cdf[idx] - BBSpecProbBelow;

  if (a != 0.)
    thisEvt.energy = x[idx] + (-b + std::sqrt(b * b - 4. * a * c)) / 2. / a;
  else
    thisEvt.energy = x[idx] - c / b;

  // ---- surface selection: weighted by effective area (area × emissivity) ----
  // Photon emission per surface scales as ε·A·T⁴ (shared T here), so a gray
  // surface must be selected — and emit — in proportion to ε·A, not A.
  if (surfaces.empty()) {
    G4Exception("ThermalSurface::GenEvt", "BBR019", FatalException,
                "No surfaces added. Call AddBoxSurface before GenEvt.");
  }
  G4double totalEffArea = 0.;
  int      N_surfaces   = (int)surfaces.size();
  for (int i = 0; i < N_surfaces; ++i)
    totalEffArea += surfaces[i].area * surfaces[i].emissivity;
  if (!(totalEffArea > 0.) || !std::isfinite(totalEffArea))
    G4Exception("ThermalSurface::GenEvt", "BBR019", FatalException,
                "No emitting area. Add a surface with positive emissivity.");

  G4double prob  = G4UniformRand() * totalEffArea;
  int idx_s = -1;
  for (int i = 0; i < N_surfaces; ++i) {
    const G4double weight = surfaces[i].area * surfaces[i].emissivity;
    if (weight == 0.) continue;
    idx_s = i;   // float residue past the end lands on the last emitting surface
    prob -= weight;
    if (prob <= 0.) break;
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
