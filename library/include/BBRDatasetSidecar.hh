#ifndef BBRDatasetSidecar_hh
#define BBRDatasetSidecar_hh

#include "G4Types.hh"

#include <string>
#include <utility>
#include <vector>

class G4VSolid;

// The HFSS dataset sidecar <waveguides>/<dirStem>.dataset.json, schema 1.x,
// agreed with Blackbody-Simulations on 2026-10-05 and described in
// validation/README.md, Dataset sidecars (field table, checks F1-F13 and
// C1-C5). Parse() reads it and runs the geometry-free checks F1-F11 and F13;
// CheckFitsSolid() is F12, run by BBRCrackLibrary::ValidatePlacedCracks; the CSV
// checks C1-C5 run in BBRHFSSData. The Python twin is bbrsim.sidecar.
//
// Frames: canonical (p, l, g) = (propagation, long, gap) = crack-local (x, y, z).
// The sampler implements one frame: HFSS global X = -g, Y = +l, Z = +p, and the
// exit CS (Z, Y, -X) = (p, l, g). A sidecar declaring another is refused (BBR025).
//
// Fatal G4Exception codes:
//   BBR024  missing, unreadable or malformed sidecar, a mistyped field, or one
//           that disagrees with its directory (ID, frequency) or, through
//           BBRCrackLibrary, with its sibling frequencies' invariant blocks
//   BBR025  a convention, frame, symmetry, cross-section or mode the sampler
//           does not implement, or a declared exit section that does not fit
//           strictly inside the placed crack solid
struct BBRCrossSection {
  enum Shape { kRectangle, kDisc };
  Shape shape = kRectangle;
  G4double yHalf_m = 0., zHalf_m = 0.;   // rectangle: half-extents along y_e (l) and z_e (g)
  G4double radius_m = 0.;                // disc
  // (y, z) [m] lies inside the section inflated by relTol: |y| <= y_half (1 + relTol)
  // (likewise z), or y^2 + z^2 <= R^2 (1 + relTol). The default 1e-6 is the HFSS
  // runner's own containment tolerance; never use a stricter one here, or a
  // producer-valid dataset is rejected.
  G4bool Contains(G4double y_m, G4double z_m, G4double relTol = 1e-6) const;
  // Points [m] on the boundary inflated by relTol: corners and edge midpoints of
  // a rectangle, nDisc points around a disc.
  std::vector<std::pair<G4double, G4double>> Boundary(G4double relTol = 1e-6, G4int nDisc = 360) const;
};

struct BBRAxisRange {
  G4double min = 0., max = 0.;
  G4int count = 0;
};

struct BBRDatasetSidecar {
  // Reads <waveguidesDir>/<dirStem>.dataset.json; dirStem = <datasetId>_<token>GHz.
  static BBRDatasetSidecar Load(const std::string& waveguidesDir, const std::string& datasetId,
                                const std::string& dirStem, G4double dirFreqGHz);
  // The same checks on JSON text; path names the file in messages.
  static BBRDatasetSidecar Parse(const std::string& text, const std::string& path,
                                 const std::string& datasetId, const std::string& dirStem,
                                 G4double dirFreqGHz);

  // F12: fatal BBR025 unless the solid's local origin is inside it and every
  // boundary point of the cross-section, at both exit faces (+-x) and under
  // both transverse mirrors, is kInside after the wrapper's axial inset.
  void CheckFitsSolid(const G4VSolid& solid, const std::string& volumeName) const;

  std::string path, datasetId, frequencyLabel;
  G4double frequencyGHz = -1.;
  std::vector<G4double> incidentPhiDeg, incidentThetaDeg;
  std::vector<std::string> farFieldColumns, exitFieldColumns;
  BBRAxisRange farFieldTheta, farFieldPhi;   // degrees
  G4int farFieldPointsPerKey = 0;
  BBRAxisRange exitY, exitZ;                 // metres
  G4int exitPointsPerKey = 0;
  std::string outsidePoints;                 // "none" | "omitted" | "zero" (C5 exempts zero-field points for "zero")
  BBRCrossSection crossSection;
  G4double extentP_mm = 0., extentL_mm = 0., extentG_mm = 0.;
  std::string lowestMode;
  G4double cutoffGHz = -1., polarizationFilterLimitGHz = -1.;
  // JSON text of the frequency-independent physics: the frame mapping
  // (hfss_global_axes_in_canonical, exit_cs_axes_in_canonical), symmetry,
  // boundaries, geometry shape and extent_mm, modes without basis (recorded
  // only) and propagating_count (per frequency), and the exit cross-section.
  // Descriptive and pose-dependent fields are left out, so sidecars of one
  // dataset from different writers or poses still agree. Compare with
  // SameInvariant, never as text.
  std::string invariant;

  // True when this and other carry the same frequency-independent physics:
  // objects compared over the keys both carry, arrays element by element,
  // numbers to 1e-9 relative, strings and booleans exactly.
  G4bool SameInvariant(const BBRDatasetSidecar& other) const;
};

#endif
