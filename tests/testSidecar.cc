// testSidecar — BBRDatasetSidecar: the schema-1 checks F1-F11 and F13 on
// sidecar text (tests/HFSSFixture.hh builds a valid one), and the fit of the
// declared exit section inside a crack solid (F12).
#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "BBRDatasetSidecar.hh"

#include "G4Box.hh"
#include "G4DisplacedSolid.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"

using namespace bbrtest;

namespace {
const std::string kStem = "crack_500GHz";
std::string Valid() { return hfssfix::SidecarJson(hfssfix::SidecarFrom(kStem, hfssfix::Mini500())); }
BBRDatasetSidecar Parse(const std::string& text, const std::string& stem = kStem, double f = 500.) {
  return BBRDatasetSidecar::Parse(text, "test.dataset.json", "crack", stem, f);
}
// text with `from` replaced by `to`; throws if `from` is absent, so a stale edit cannot pass.
std::string Edit(std::string text, const std::string& from, const std::string& to) {
  const auto p = text.find(from);
  if (p == std::string::npos) throw std::runtime_error("Edit: '" + from + "' not in the sidecar text");
  return text.replace(p, from.size(), to);
}
// Parse(text) must stop with `code` and a message containing `what`, so each case
// pins the check it is meant to reach, not merely one that fires earlier.
void ExpectFail(const std::string& code, const std::string& what, const std::string& text) {
  try {
    Parse(text);
    Report(false, __FILE__, __LINE__, "expected " + code + " (" + what + "), none raised");
  } catch (const G4ExceptionCaught& e) {
    Report(e.code == code && e.origin == "BBRDatasetSidecar" && std::string(e.what()).find(what) != std::string::npos,
           __FILE__, __LINE__, "expected " + code + " containing \"" + what + "\", got " + e.code + " from " + e.origin + ": " + e.what());
  } catch (const std::exception& e) {
    Report(false, __FILE__, __LINE__, "expected " + code + ", got C++ exception: " + e.what());
  }
}
G4VSolid* Tube(double radius) {   // a G4Tubs turned so its axis is local x, as the wrapper needs
  G4RotationMatrix ry;
  ry.rotateY(90. * deg);
  return new G4DisplacedSolid("tube", new G4Tubs("t", 0., radius, 0.2 * mm, 0., CLHEP::twopi),
                              G4Transform3D(ry, G4ThreeVector()));
}
}  // namespace

int main(int argc, char** argv) {
  return RunCase(argc, argv, {
    {"parses_fixture", [] {
      const auto sc = Parse(Valid());
      CHECK(sc.datasetId == "crack" && sc.frequencyLabel == "500GHz");
      CHECK_NEAR(sc.frequencyGHz, 500., 0);
      CHECK(sc.crossSection.shape == BBRCrossSection::kRectangle);
      CHECK_REL(sc.crossSection.yHalf_m, 4.5e-3, 1e-12);
      CHECK(sc.lowestMode == "TE10");
      CHECK_REL(sc.cutoffGHz, 299792458. / (2 * 9e-3) / 1e9, 1e-12);
      CHECK(sc.farFieldPointsPerKey == 1 && sc.exitPointsPerKey == 1);
      CHECK(sc.incidentThetaDeg.size() == 1 && sc.incidentThetaDeg[0] == 180.);
      CHECK(!sc.invariant.empty());
      // Load reads <dir>/<stem>.dataset.json and names it by its absolute path, also from a relative dir.
      TempDir td;
      WriteFile(td.path() / "waveguides" / (kStem + ".dataset.json"), Valid());
      const auto la = BBRDatasetSidecar::Load((td.path() / "waveguides").string(), "crack", kStem, 500.);
      CHECK(la.lowestMode == "TE10" && la.SameInvariant(sc));
      CHECK(std::filesystem::path(la.path).is_absolute());
      CHECK(std::filesystem::equivalent(la.path, td.path() / "waveguides" / (kStem + ".dataset.json")));
      const auto cwd = std::filesystem::current_path();
      std::filesystem::current_path(td.path());
      const auto lr = BBRDatasetSidecar::Load("waveguides", "crack", kStem, 500.);
      std::filesystem::current_path(cwd);
      CHECK(std::filesystem::path(lr.path).is_absolute());
      CHECK(std::filesystem::equivalent(lr.path, td.path() / "waveguides" / (kStem + ".dataset.json")));
    }},
    {"forward_compatible", [] {
      // A later 1.x with unknown fields, written with CRLF line ends and a UTF-8 BOM, parses.
      const std::string t = Edit(Valid(), "\"schema_version\": \"1.0\"",
                                 "\"schema_version\": \"1.3\", \"x_future\": {\"a\": [1, 2]}");
      std::string crlf;
      for (char ch : t) crlf += (ch == '\n') ? std::string("\r\n") : std::string(1, ch);
      CHECK(Parse("\xEF\xBB\xBF" + crlf).lowestMode == "TE10");
    }},
    {"f1_f3_identity", [] {
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"1.0\"", "\"2.0\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"dataset_id\": \"crack\"", "\"dataset_id\": \"other\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { Parse(Valid(), "crack_501GHz", 501.); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"frequency_ghz\": 500", "\"frequency_ghz\": 520")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": \"1\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { Parse("{not json"); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR024", [] { BBRDatasetSidecar::Load("/nonexistent", "crack", kStem, 500.); }, "BBRDatasetSidecar");
    }},
    {"f4_frames", [] {
      const std::string hfssAxes = "{\"x\": [0, 0, -1], \"y\": [0, 1, 0], \"z\": [1, 0, 0]}";
      const std::string exitAxes = "\"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, 1, 0], \"z\": [0, 0, 1]}";
      // A left-handed HFSS axis map.
      ExpectFail("BBR025", "not a right-handed orthonormal basis",
                 Edit(Valid(), hfssAxes, "{\"x\": [0, 0, 1], \"y\": [0, 1, 0], \"z\": [1, 0, 0]}"));
      // Inconsistent exit_cs mapping: the HFSS axes rotated 180 degrees about p, exit_cs_axes_in_canonical left as is.
      const std::string rotated = Edit(Valid(), hfssAxes, "{\"x\": [0, 0, 1], \"y\": [0, -1, 0], \"z\": [1, 0, 0]}");
      ExpectFail("BBR025", "exit_cs_axes_in_canonical disagrees", rotated);
      // The same rotation labelled consistently (the registry's first version): a proper,
      // self-consistent frame, but not the one the sampler implements.
      ExpectFail("BBR025", "the frames differ from the one the sampler implements",
                 Edit(rotated, exitAxes, "\"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, -1, 0], \"z\": [0, 0, -1]}"));
      ExpectFail("BBR025", "is not x cross y", Edit(Valid(), "\"z\": [-1, 0, 0]}", "\"z\": [1, 0, 0]}"));
      ExpectFail("BBR025", "outward normal must map to +p",
                 Edit(Valid(), "\"exit_outward_normal_global\": [0, 0, 1]", "\"exit_outward_normal_global\": [0, 0, -1]"));
    }},
    {"f5_f10_conventions", [] {
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"arrival_direction\"", "\"propagation_direction\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"field_components_frame\": \"hfss_global\"", "\"field_components_frame\": \"exit_cs\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"shape\": \"rectangle\"", "\"shape\": \"polygon\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"incoming_includes_cos_theta\": false", "\"incoming_includes_cos_theta\": true")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"mirror_l\": true", "\"mirror_l\": false")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"Phi\", \"Theta\", \"rEphi_real\"", "\"Theta\", \"Phi\", \"rEphi_real\"")); }, "BBRDatasetSidecar");
    }},
    {"f13_modes", [] {
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.yHalf = 5e-3; p.zHalf = 2.5e-5;                 // crack1's HFSS section
      const auto box = Parse(hfssfix::SidecarJson(p));
      CHECK_NEAR(box.cutoffGHz, 14.9896229, 1e-6);
      CHECK_NEAR(box.polarizationFilterLimitGHz, 2997.92458, 1e-4);
      p.disc = true; p.radius = 5e-5;                   // the round gap
      const auto disc = Parse(hfssfix::SidecarJson(p));
      CHECK(disc.lowestMode == "TE11");
      CHECK_NEAR(disc.cutoffGHz, 1756.98, 0.01);
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"mode\": \"TE10\"", "\"mode\": \"TE01\"")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"cutoff_ghz\": ", "\"cutoff_ghz\": 1")); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"polarization_filter_limit_ghz\": ",
                                                  "\"polarization_filter_limit_ghz\": 1")); }, "BBRDatasetSidecar");
      // A section taller along g than along l: the lowest mode is TE01 at c/(2b), which is also the filter limit.
      auto q = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      q.yHalf = 1e-5; q.zHalf = 2e-5;
      const auto tall = Parse(hfssfix::SidecarJson(q));
      CHECK(tall.lowestMode == "TE01");
      CHECK_REL(tall.cutoffGHz, 299792458. / (2 * 4e-5) / 1e9, 1e-12);
      CHECK_REL(tall.polarizationFilterLimitGHz, 299792458. / (2 * 4e-5) / 1e9, 1e-12);
    }},
    {"invariant_physics_only", [] {
      // Descriptive and pose fields differ between writers: same physics.
      const auto a = Parse(Valid());
      std::string t = Edit(Valid(), "\"canonical\": \"p,l,g\"", "\"canonical\": \"another writer's wording\"");
      t = Edit(t, "\"pose_rule\": \"canonical-z\"", "\"pose_rule\": \"legacy\"");
      const auto b = Parse(t);
      CHECK(a.SameInvariant(b) && b.SameInvariant(a));
      // propagating_count (per frequency) and basis (recorded only) are not compared.
      const std::string c33 = Edit(Valid(), "\"propagating_count\": 0", "\"propagating_count\": 33");
      const auto c = Parse(Edit(c33, "\"modes\": {", "\"modes\": {\"basis\": \"closed PEC a=0.009 b=5e-05\", "));
      const auto d = Parse(Edit(Valid(), "\"modes\": {", "\"modes\": {\"basis\": \"another basis string\", "));
      CHECK(a.SameInvariant(c) && c.SameInvariant(a));
      CHECK(c.SameInvariant(d) && d.SameInvariant(c));
      // A sidecar not produced by Parse (empty invariant) matches nothing and does not throw.
      const BBRDatasetSidecar blank;
      CHECK(!blank.SameInvariant(a) && !a.SameInvariant(blank) && !blank.SameInvariant(blank));
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.zHalf = 2.6e-5;                                   // another exit section: different physics
      CHECK(!a.SameInvariant(Parse(hfssfix::SidecarJson(p))));
    }},
    {"field_types", [] {
      // Mistyped fields are BBR024; a column list that is an array but not the
      // positional order (non-string entries included) is BBR025, as in bbrsim.sidecar.check.
      ExpectFail("BBR024", "boundaries must be an object",
                 Edit(Valid(), "\"boundaries\": {\"entrance\": \"radiation\", \"exit\": \"radiation\", \"walls\": \"PEC\"}",
                      "\"boundaries\": \"PEC\""));
      ExpectFail("BBR024", "far_field.points_per_key must be an integer",
                 Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": 4294967297"));
      ExpectFail("BBR024", "far_field.points_per_key must be an integer",
                 Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": -4294967297"));
      ExpectFail("BBR024", "far_field.points_per_key must be an integer",
                 Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": 1.5"));
      ExpectFail("BBR024", "far_field.theta_deg.count must be an integer",
                 Edit(Valid(), "\"count\": 1", "\"count\": 2147483648"));
      ExpectFail("BBR024", "far_field.theta_deg needs min <= max",
                 Edit(Valid(), "\"count\": 1", "\"count\": 0"));
      const std::string ffCols = "\"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"Phi\"";
      const std::string xfCols = "\"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"OutgoingPower\"";
      ExpectFail("BBR024", "far_field.columns must be an array", Edit(Valid(), ffCols, "\"columns\": \"Freq\", \"x\": [\"Phi\""));
      ExpectFail("BBR024", "exit_field.columns must be an array", Edit(Valid(), xfCols, "\"columns\": {}, \"x\": [\"OutgoingPower\""));
      ExpectFail("BBR025", "far_field.columns differ from the positional order",
                 Edit(Valid(), ffCols, "\"columns\": [1, \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"Phi\""));
      ExpectFail("BBR025", "exit_field.columns differ from the positional order",
                 Edit(Valid(), xfCols, "\"columns\": [\"Freq\", null, \"IWavePhi\", \"IWaveTheta\", \"OutgoingPower\""));
      ExpectFail("BBR025", "exit_field.columns differ from the positional order",
                 Edit(Valid(), xfCols, "\"columns\": [\"Freq\", \"IWavePhi\", \"Ephi\", \"IWaveTheta\", \"OutgoingPower\""));
    }},
    {"f12_fits_solid", [] {
      const auto box = Parse(Valid());                  // section 4.5 mm x 25 um
      box.CheckFitsSolid(G4Box("ok", 2 * mm, 5 * mm, 0.026 * mm), "ok");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("thin", 2 * mm, 5 * mm, 0.025 * mm), "thin"); }, "BBRDatasetSidecar");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("short", 2 * mm, 4 * mm, 0.026 * mm), "short"); }, "BBRDatasetSidecar");
      // A solid whose local origin is not inside it: the wrapper finds the exit face from the origin.
      G4DisplacedSolid off("off", new G4Box("b", 2 * mm, 5 * mm, 0.026 * mm), G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 1 * mm)));
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(off, "off"); }, "BBRDatasetSidecar");
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5;
      const auto disc = Parse(hfssfix::SidecarJson(p));
      disc.CheckFitsSolid(*Tube(0.051 * mm), "tube51");
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(*Tube(0.050 * mm), "tube50"); }, "BBRDatasetSidecar");   // rim on the wall
      // A plain G4Tubs keeps its axis on local z, so local x is radial: rejected at startup.
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(G4Tubs("bare", 0., 0.051 * mm, 0.2 * mm, 0., CLHEP::twopi), "bare"); }, "BBRDatasetSidecar");
    }},
  });
}
