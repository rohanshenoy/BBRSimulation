#ifndef BBRDatasetSidecar_hh
#define BBRDatasetSidecar_hh

#include "G4Types.hh"

#include <string>
#include <utility>
#include <vector>

class G4VSolid;

// The HFSS dataset sidecar <waveguides>/<dirStem>.dataset.json, schema 1.x or
// 2.x (OR-034, OR-036, HF-047), agreed with Blackbody-Simulations (1.x on
// 2026-10-05) and described in
// validation/README.md, Dataset sidecars (field table, checks F1-F13 and
// C1-C5). Parse() reads it and runs the geometry-free checks F1-F11 and F13
// (F11 includes a rectangle's extent_mm l and g against its declared section);
// CheckFitsSolid() is F12, run by BBRCrackLibrary::ValidatePlacedCracks: the
// declared section and extent_mm.p against the placed solid; the CSV
// checks C1-C5 run in BBRHFSSData. The Python twin is bbrsim.sidecar.
//
// Frames: canonical (p, l, g) = (propagation, long, gap) = crack-local (x, y, z).
// The sampler implements one frame: HFSS global X = -g, Y = +l, Z = +p, and the
// exit CS (Z, Y, -X) = (p, l, g). A sidecar declaring another is refused (BBR025).
//
// Fatal G4Exception codes:
//   BBR024  missing, unreadable or malformed sidecar, a mistyped field, or one
//           that disagrees with its directory (ID, frequency) or, through
//           BBRCrackLibrary, with its sibling frequencies' invariant blocks;
//           under 2.x also a per-key list that does not match the declared
//           keys exactly, a lattice step that is not (max - min)/(count - 1),
//           and a dark basis whose t_flux is not exactly 0 (OR-038)
//   BBR025  a convention, frame, symmetry, cross-section or mode the sampler
//           does not implement, or a declared exit section that does not equal
//           the placed crack solid's section (F12; OR-020); under 2.x also a
//           fixed definition string that differs, a c_k or cell_90 that is not
//           the closed form of the declared grid, and an emitter assumption
//           other than uniform theta
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
  // a rectangle, nDisc points around a disc. offset_m is added to each half-size, or
  // to the radius, after the relTol inflation; negative moves the points inward.
  std::vector<std::pair<G4double, G4double>> Boundary(G4double relTol = 1e-6, G4int nDisc = 360,
                                                      G4double offset_m = 0.) const;
  static constexpr G4double kClampInset_m = 1e-9;   // OR-020: exit points are clamped this far inside the section
};

struct BBRAxisRange {
  G4double min = 0., max = 0.;
  G4int count = 0;
  G4double step = 0.;   // 2.x only: the declared lattice step, equal to (max - min)/(count - 1); 0 under 1.x
};

// One declared incidence key of a 2.x sidecar, matched across its per-key lists
// (transmittance.c_key and per_key, exit_field.scale.per_key, far_field.ff_over_flux).
struct BBRKeyRecord {
  G4double phiDeg = 0., thetaDeg = 0.;
  G4double cK = 0.;                // the key's InterceptFactor divisor
  G4double tFlux0 = 0., tFlux1 = 0.;   // stored OutgoingPower / (A I) per polarization, before the c_k division
  G4double w0 = 1., w1 = 1.;       // exit_field.scale per-key factors
  G4bool dark0 = false, dark1 = false;   // a dark basis: T exactly 0, never sampled (OR-038)
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

  // HF-049: the divisor of the theta = 90 deg key, 1/<1/max(cos a, 0.01)> for crossings uniform in a
  // over the declared grid's outermost band, a from 90 - Delta/2 to 90 deg, Delta = (the smallest declared
  // theta above 90) - 90. -1 when the list has no 90 deg entry (within 1e-9) or nothing above it.
  static G4double Cell90(const std::vector<G4double>& thetaDeg);
  // c_k of a key (OR-034): |cos theta|, or Cell90(thetaList) at theta = 90 deg.
  static G4double InterceptFactor(G4double thetaDeg, const std::vector<G4double>& thetaList);
  // The 2.x record of a declared key, matched within 1e-9 in both angles; nullptr if absent (and
  // always under 1.x, which has no per-key records).
  const BBRKeyRecord* Key(G4double phiDeg, G4double thetaDeg) const;

  // The regime-changing cutoffs of a section, ascending (OR-031): a rectangle's TE10 at c/(2a), a = 2 y_e_half,
  // then TE0n at n c/(2b), b = 2 z_e_half, for n = 1, 2, ... up to 20000 GHz inclusive, labelled TE{m}{n} when
  // both indices are below 10 and TE{m},{n} otherwise; a disc's TE11, TM01 and TE21. Closed forms, no Bessel
  // function: the zeros x'11, x01 and x'21 are constants.
  static std::vector<std::pair<std::string, G4double>> RegimeCutoffs(const BBRCrossSection& section);

  G4int schemaMajor = 1;               // 1 or 2, from schema_version
  std::vector<BBRKeyRecord> keys;      // 2.x only: one per declared key, phi-major then theta
  G4double cell90 = -1.;               // 2.x only: transmittance.cell_90

  std::string path, datasetId, frequencyLabel;
  G4double frequencyGHz = -1.;
  std::vector<G4double> incidentPhiDeg, incidentThetaDeg;
  std::vector<std::string> farFieldColumns, exitFieldColumns;   // exitFieldColumns: the 1.0 or the 2.0 list
  BBRAxisRange farFieldTheta, farFieldPhi;   // degrees
  G4int farFieldPointsPerKey = 0;
  BBRAxisRange exitY, exitZ;                 // metres
  G4int exitPointsPerKey = 0;
  std::string outsidePoints;                 // "none" | "omitted" | "zero" (C5 exempts zero-field points for "zero")
  BBRCrossSection crossSection;
  G4double extentP_mm = 0., extentL_mm = 0., extentG_mm = 0.;
  std::string lowestMode;
  G4double cutoffGHz = -1., polarizationFilterLimitGHz = -1.;
  // Every cutoff a photon can cross to reach a table of this dataset (BBR026): 2.x the declared
  // modes.regime_cutoffs_ghz, checked against RegimeCutoffs of the section; 1.x {{lowestMode, cutoffGHz}}.
  std::vector<std::pair<std::string, G4double>> regimeCutoffs;
  G4double rhoValidBelowGHz = -1.;     // 2.x only: cross_term.valid_below_ghz, equal to polarizationFilterLimitGHz
  // JSON text of the frequency-independent physics: the frame mapping
  // (hfss_global_axes_in_canonical, exit_cs_axes_in_canonical), symmetry,
  // boundaries, geometry shape and extent_mm, modes without basis (recorded
  // only) and propagating_count (per frequency), and the exit cross-section.
  // Descriptive and pose-dependent fields are left out, so sidecars of one
  // dataset from different writers or poses still agree. Compare with
  // SameInvariant, never as text.
  std::string invariant;
  // Compact JSON text of the recorded-only fields, copied verbatim and never
  // checked, for <stem>.metadata.json (BBRCrackLibrary::PlacedCracksJson):
  // {"provenance", "boundaries", "geometry": {"shape", "bounding_box_mm"},
  //  "symmetry": {"rotational"}, "frames": {"pose_rule", "exit_cs_origin_mm_global"
  //  (frames.exit_cs.origin_mm_global)}, "excitation": {"origin_mm_global"}}.
  // A field the sidecar lacks is null.
  std::string recorded;

  // The top-level invariant blocks on which this and other differ, in the order
  // schema, excitation, frames, symmetry, boundaries, geometry, modes, cross_section (schema is the major
  // version, excitation the incidence grid: one ID must not mix either across frequencies, G4-032): objects
  // compared over the keys both carry, arrays element by element, numbers to
  // 1e-9 relative, strings and booleans exactly. Empty when either invariant
  // is empty (a sidecar not produced by Parse). The Python twin is
  // bbrsim.sidecar.invariant_diff.
  std::vector<std::string> InvariantDiff(const BBRDatasetSidecar& other) const;

  // True when this and other carry the same frequency-independent physics:
  // both produced by Parse and InvariantDiff(other) empty.
  G4bool SameInvariant(const BBRDatasetSidecar& other) const;
};

#endif
