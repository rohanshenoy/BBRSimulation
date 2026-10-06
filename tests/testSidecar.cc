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

#include <nlohmann/json.hpp>

#include <map>

using namespace bbrtest;
using json = nlohmann::json;

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
      // Each message substring pins the check the line names, so no case can pass on an earlier
      // check that raises the same code.
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"1.0\"", "\"2.0\"")); }, "BBRDatasetSidecar",
                        "schema_version 2.0: BBRsim reads 1.x");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"dataset_id\": \"crack\"", "\"dataset_id\": \"other\"")); }, "BBRDatasetSidecar",
                        "dataset_id \"other\", but the directory and crack volume say \"crack\"");
      ExpectG4Exception("BBR024", [] { Parse(Valid(), "crack_501GHz", 501.); }, "BBRDatasetSidecar",
                        "frequency_label \"500GHz\", but the directory says \"501GHz\"");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"frequency_ghz\": 500", "\"frequency_ghz\": 520")); }, "BBRDatasetSidecar",
                        "frequency_ghz 520 disagrees with the directory frequency 500 GHz (0.1 %)");
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": \"1\"")); }, "BBRDatasetSidecar",
                        "far_field.points_per_key must be an integer within the int range");
      ExpectG4Exception("BBR024", [] { Parse("{not json"); }, "BBRDatasetSidecar", "not valid JSON: ");
      ExpectG4Exception("BBR024", [] { BBRDatasetSidecar::Load("/nonexistent", "crack", kStem, 500.); }, "BBRDatasetSidecar",
                        "cannot be opened. Every HFSS dataset needs this sidecar");
    }},
    {"f4_frames", [] {
      const std::string hfssAxes = "{\"x\": [0, 0, -1], \"y\": [0, 1, 0], \"z\": [1, 0, 0]}";
      const std::string exitAxes = "\"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, 1, 0], \"z\": [0, 0, 1]}";
      // A left-handed HFSS axis map.
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), hfssAxes, "{\"x\": [0, 0, 1], \"y\": [0, 1, 0], \"z\": [1, 0, 0]}")); },
                        "BBRDatasetSidecar", "not a right-handed orthonormal basis");
      // Inconsistent exit_cs mapping: the HFSS axes rotated 180 degrees about p, exit_cs_axes_in_canonical left as is.
      const std::string rotated = Edit(Valid(), hfssAxes, "{\"x\": [0, 0, 1], \"y\": [0, -1, 0], \"z\": [1, 0, 0]}");
      ExpectG4Exception("BBR025", [&] { Parse(rotated); }, "BBRDatasetSidecar", "exit_cs_axes_in_canonical disagrees");
      // The same rotation labelled consistently (the registry's first version): a proper,
      // self-consistent frame, but not the one the sampler implements.
      ExpectG4Exception("BBR025", [&] { Parse(Edit(rotated, exitAxes, "\"exit_cs_axes_in_canonical\": {\"x\": [1, 0, 0], \"y\": [0, -1, 0], \"z\": [0, 0, -1]}")); },
                        "BBRDatasetSidecar", "the frames differ from the one the sampler implements");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), "\"z\": [-1, 0, 0]}", "\"z\": [1, 0, 0]}")); },
                        "BBRDatasetSidecar", "is not x cross y");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), "\"exit_outward_normal_global\": [0, 0, 1]", "\"exit_outward_normal_global\": [0, 0, -1]")); },
                        "BBRDatasetSidecar", "outward normal must map to +p");
    }},
    {"f5_f10_conventions", [] {
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"arrival_direction\"", "\"propagation_direction\"")); }, "BBRDatasetSidecar",
                        "excitation.incidence_convention is \"propagation_direction\"; BBRsim implements only \"arrival_direction\"");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"field_components_frame\": \"hfss_global\"", "\"field_components_frame\": \"exit_cs\"")); }, "BBRDatasetSidecar",
                        "exit_field.field_in_ref_cs and field_components_frame contradict each other");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"shape\": \"rectangle\"", "\"shape\": \"polygon\"")); }, "BBRDatasetSidecar",
                        "cross_section shape polygon is reserved; BBRsim does not support it yet");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"incoming_includes_cos_theta\": false", "\"incoming_includes_cos_theta\": true")); }, "BBRDatasetSidecar",
                        "transmittance.incoming_includes_cos_theta must be false");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"mirror_l\": true", "\"mirror_l\": false")); }, "BBRDatasetSidecar",
                        "symmetry.mirror_l is false; the sampler folds by both transverse mirrors");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"Phi\", \"Theta\", \"rEphi_real\"", "\"Theta\", \"Phi\", \"rEphi_real\"")); }, "BBRDatasetSidecar",
                        "far_field.columns differ from the positional order BBRHFSSData reads");
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
      // The message prints the declared triple, then the re-derived one, so the substring names
      // the field that was edited (the fixture's 9 mm x 50 um section: TE10, 16.655..., 2997.92458 GHz).
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"mode\": \"TE10\"", "\"mode\": \"TE01\"")); }, "BBRDatasetSidecar",
                        "modes (TE01, 16.6551365556, 2997.92458 GHz) disagree with the cross-section (TE10, 16.6551365556, 2997.92458 GHz)");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"cutoff_ghz\": ", "\"cutoff_ghz\": 1")); }, "BBRDatasetSidecar",
                        "modes (TE10, 116.655136556, 2997.92458 GHz) disagree with the cross-section (TE10, 16.6551365556, 2997.92458 GHz)");
      ExpectG4Exception("BBR025", [] { Parse(Edit(Valid(), "\"polarization_filter_limit_ghz\": ",
                                                  "\"polarization_filter_limit_ghz\": 1")); }, "BBRDatasetSidecar",
                        "modes (TE10, 16.6551365556, 12997.92458 GHz) disagree with the cross-section (TE10, 16.6551365556, 2997.92458 GHz)");
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
    {"invariant_diff", [] {
      // The names of the differing top-level blocks, in the order frames, symmetry,
      // boundaries, geometry, modes, cross_section; empty exactly when SameInvariant.
      using V = std::vector<std::string>;
      const auto a = Parse(Valid());
      const auto pose = Parse(Edit(Valid(), "\"pose_rule\": \"canonical-z\"", "\"pose_rule\": \"legacy\""));
      CHECK(a.InvariantDiff(pose).empty() && a.SameInvariant(pose));
      const auto sym = Parse(Edit(Valid(), "\"rotational\": false", "\"rotational\": true"));
      CHECK(a.InvariantDiff(sym) == V{"symmetry"} && !a.SameInvariant(sym));
      const auto wall = Parse(Edit(Valid(), "\"walls\": \"PEC\"", "\"walls\": \"Cu\""));
      CHECK(a.InvariantDiff(wall) == V{"boundaries"});
      const auto both = Parse(Edit(Edit(Valid(), "\"walls\": \"PEC\"", "\"walls\": \"Cu\""),
                                   "\"rotational\": false", "\"rotational\": true"));
      CHECK(both.InvariantDiff(a) == (V{"symmetry", "boundaries"}));
      const auto geo = Parse(Edit(Valid(), "\"g\": 0.052", "\"g\": 0.06"));
      CHECK(a.InvariantDiff(geo) == V{"geometry"});
      const auto near = Parse(Edit(Valid(), "\"g\": 0.052", "\"g\": 0.0520000000001"));   // 2e-12 relative
      CHECK(a.InvariantDiff(near).empty());
      // Another exit section moves the modes it determines too (TE10 follows the long side).
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.yHalf = 4.0e-3;
      const auto sec = Parse(hfssfix::SidecarJson(p));
      CHECK(a.InvariantDiff(sec) == (V{"modes", "cross_section"}) && sec.InvariantDiff(a) == a.InvariantDiff(sec));
      // A sidecar not produced by Parse has no invariant: no blocks named, and SameInvariant stays false.
      const BBRDatasetSidecar blank;
      CHECK(blank.InvariantDiff(a).empty() && a.InvariantDiff(blank).empty() && !a.SameInvariant(blank));
    }},
    {"field_types", [] {
      // Mistyped fields are BBR024; a column list that is an array but not the
      // positional order (non-string entries included) is BBR025, as in bbrsim.sidecar.check.
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"boundaries\": {\"entrance\": \"radiation\", \"exit\": \"radiation\", \"walls\": \"PEC\"}",
                                                   "\"boundaries\": \"PEC\"")); },
                        "BBRDatasetSidecar", "boundaries must be an object");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": 4294967297")); },
                        "BBRDatasetSidecar", "far_field.points_per_key must be an integer");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": -4294967297")); },
                        "BBRDatasetSidecar", "far_field.points_per_key must be an integer");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"points_per_key\": 1", "\"points_per_key\": 1.5")); },
                        "BBRDatasetSidecar", "far_field.points_per_key must be an integer");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"count\": 1", "\"count\": 2147483648")); },
                        "BBRDatasetSidecar", "far_field.theta_deg.count must be an integer");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), "\"count\": 1", "\"count\": 0")); },
                        "BBRDatasetSidecar", "far_field.theta_deg needs min <= max");
      const std::string ffCols = "\"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"Phi\"";
      const std::string xfCols = "\"columns\": [\"Freq\", \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"OutgoingPower\"";
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), ffCols, "\"columns\": \"Freq\", \"x\": [\"Phi\"")); },
                        "BBRDatasetSidecar", "far_field.columns must be an array");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), xfCols, "\"columns\": {}, \"x\": [\"OutgoingPower\"")); },
                        "BBRDatasetSidecar", "exit_field.columns must be an array");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), ffCols, "\"columns\": [1, \"Ephi\", \"IWavePhi\", \"IWaveTheta\", \"Phi\"")); },
                        "BBRDatasetSidecar", "far_field.columns differ from the positional order");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), xfCols, "\"columns\": [\"Freq\", null, \"IWavePhi\", \"IWaveTheta\", \"OutgoingPower\"")); },
                        "BBRDatasetSidecar", "exit_field.columns differ from the positional order");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), xfCols, "\"columns\": [\"Freq\", \"IWavePhi\", \"Ephi\", \"IWaveTheta\", \"OutgoingPower\"")); },
                        "BBRDatasetSidecar", "exit_field.columns differ from the positional order");
    }},
    {"recorded_fields", [] {
      // The recorded-only fields, copied verbatim for metadata.json; absent ones are null.
      const auto r = json::parse(Parse(Valid()).recorded);
      CHECK(r.at("frames").at("pose_rule") == "canonical-z");
      CHECK(r.at("geometry").at("shape") == "box");
      CHECK(r.at("symmetry").at("rotational") == false);
      CHECK(r.at("excitation").at("origin_mm_global").is_array() && r.at("excitation").at("origin_mm_global").size() == 3);
      CHECK(r.at("frames").at("exit_cs_origin_mm_global").is_array() && r.at("frames").at("exit_cs_origin_mm_global").size() == 3);
      CHECK(r.at("geometry").at("bounding_box_mm").is_array() && r.at("geometry").at("bounding_box_mm").size() == 6);
      const json source = json::parse(Valid());
      CHECK(r.at("provenance") == source.at("provenance") && r.at("boundaries") == source.at("boundaries"));
      const std::map<std::string, std::vector<std::string>> keys = {
          {"provenance", {}}, {"boundaries", {}}, {"geometry", {"shape", "bounding_box_mm"}},
          {"symmetry", {"rotational"}}, {"frames", {"pose_rule", "exit_cs_origin_mm_global"}},
          {"excitation", {"origin_mm_global"}}};
      CHECK(r.is_object() && r.size() == keys.size());
      for (const auto& [block, fields] : keys) {
        CHECK(r.contains(block));
        if (!fields.empty()) CHECK(r.at(block).is_object() && r.at(block).size() == fields.size());
        for (const auto& f : fields) CHECK(r.at(block).contains(f));
      }
      // Fields the sidecar lacks become null; an unknown provenance field is copied verbatim.
      std::string t = Edit(Valid(), ", \"pose_rule\": \"canonical-z\"", "");
      t = Edit(t, ",\n    \"bounding_box_mm\": [-0.026, -5, 0, 0.026, 5, 4]", "");
      t = Edit(t, "\"provenance\": {", "\"provenance\": {\"x_future\": {\"a\": [1, 2]}, ");
      const auto m = json::parse(Parse(t).recorded);
      CHECK(m.at("frames").at("pose_rule").is_null());
      CHECK(m.at("geometry").at("bounding_box_mm").is_null());
      CHECK(m.at("geometry").at("shape") == "box");
      CHECK(m.at("provenance").at("x_future") == json::parse("{\"a\": [1, 2]}"));
      // An absent parent object gives null too, and Parse does not fail on it.
      const auto n = json::parse(Parse(Edit(Valid(), "\"provenance\": {\"producer\": \"tests/HFSSFixture.hh\"},\n", "")).recorded);
      CHECK(n.at("provenance").is_null() && n.at("frames").at("pose_rule") == "canonical-z");
    }},
    {"f12_fits_solid", [] {
      const auto box = Parse(Valid());                  // section 4.5 mm x 25 um
      box.CheckFitsSolid(G4Box("ok", 2 * mm, 5 * mm, 0.026 * mm), "ok");
      // The fit message names the volume; the origin message is the other check of CheckFitsSolid.
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("thin", 2 * mm, 5 * mm, 0.025 * mm), "thin"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not fit strictly inside crack volume thin: ");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("short", 2 * mm, 4 * mm, 0.026 * mm), "short"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not fit strictly inside crack volume short: ");
      // A solid whose local origin is not inside it: the wrapper finds the exit face from the origin.
      G4DisplacedSolid off("off", new G4Box("b", 2 * mm, 5 * mm, 0.026 * mm), G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 1 * mm)));
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(off, "off"); }, "BBRDatasetSidecar",
                        "the local origin of crack volume off is not inside its solid");
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5;
      const auto disc = Parse(hfssfix::SidecarJson(p));
      disc.CheckFitsSolid(*Tube(0.051 * mm), "tube51");
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(*Tube(0.050 * mm), "tube50"); }, "BBRDatasetSidecar",   // rim on the wall
                        "the declared exit cross-section does not fit strictly inside crack volume tube50: ");
      // A plain G4Tubs keeps its axis on local z, so local x is radial: rejected at startup.
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(G4Tubs("bare", 0., 0.051 * mm, 0.2 * mm, 0., CLHEP::twopi), "bare"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not fit strictly inside crack volume bare: ");
    }},
  });
}
