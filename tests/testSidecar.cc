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
// Valid() with geometry.extent_mm.p set to `p_mm` (p is compared with the solid in F12, not at parse).
std::string WithP(double p_mm, std::string text = Valid()) {
  json j = json::parse(text);
  j["geometry"]["extent_mm"]["p"] = p_mm;
  return j.dump(2);
}
// --- schema 2.0 ---
std::string Valid2() { return hfssfix::SidecarJson2(hfssfix::NineKeyParams2(kStem)); }
// The 2.0 fixture with `fn` applied to its JSON.
template <class F> std::string Mut2(F&& fn) {
  json j = json::parse(Valid2());
  fn(j);
  return j.dump(2);
}
json& At(json& j, const char* pointer) { return j[json::json_pointer(pointer)]; }
// text with the sentinel replaced, to write what nlohmann cannot dump (1e999 reads back as inf).
std::string Sentinel(std::string text, const std::string& to) {
  const std::string from = "12345.6789";
  const auto p = text.find(from);
  if (p == std::string::npos) throw std::runtime_error("Sentinel absent");
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
      ExpectG4Exception("BBR024", [] { Parse(Edit(Valid(), "\"1.0\"", "\"3.0\"")); }, "BBRDatasetSidecar",
                        "schema_version 3.0: BBRsim reads 1.x and 2.x");
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
      // boundaries can no longer differ between sidecars that parse: F10 accepts one model (CDX-003).
      const auto geo = Parse(WithP(5.));
      CHECK(a.InvariantDiff(geo) == V{"geometry"});
      const auto both = Parse(Edit(WithP(5.), "\"rotational\": false", "\"rotational\": true"));
      CHECK(both.InvariantDiff(a) == (V{"symmetry", "geometry"}));
      const auto near = Parse(WithP(4. * (1 + 2e-12)));   // 2e-12 relative
      CHECK(a.InvariantDiff(near).empty());
      // Another exit section moves the modes it determines too (TE10 follows the long side), and the
      // extent_mm l that the parse ties to it.
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.yHalf = 4.0e-3;
      const auto sec = Parse(hfssfix::SidecarJson(p));
      CHECK(a.InvariantDiff(sec) == (V{"geometry", "modes", "cross_section"}) && sec.InvariantDiff(a) == a.InvariantDiff(sec));
      // A sidecar not produced by Parse has no invariant: no blocks named, and SameInvariant stays false.
      const BBRDatasetSidecar blank;
      CHECK(blank.InvariantDiff(a).empty() && a.InvariantDiff(blank).empty() && !a.SameInvariant(blank));
    }},
    {"f10_boundaries", [] {
      // CDX-003: radiation end faces and PEC walls, the model the transport implements.
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), "\"walls\": \"PEC\"", "\"walls\": \"finite-conductivity\"")); },
                        "BBRDatasetSidecar", "boundaries.walls is \"finite-conductivity\"; BBRsim implements only \"PEC\"");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), "\"entrance\": \"radiation\"", "\"entrance\": \"port\"")); },
                        "BBRDatasetSidecar", "boundaries.entrance is \"port\"; BBRsim implements only \"radiation\"");
      ExpectG4Exception("BBR025", [&] { Parse(Edit(Valid(), "\"exit\": \"radiation\"", "\"exit\": \"PML\"")); },
                        "BBRDatasetSidecar", "boundaries.exit is \"PML\"");
      ExpectG4Exception("BBR024", [&] { Parse(Edit(Valid(), ", \"walls\": \"PEC\"", "")); },
                        "BBRDatasetSidecar", "missing field boundaries.walls");
      Parse(Valid());   // the supported declaration passes
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
    {"s2_parses_fixture", [] {
      const auto sc = Parse(Valid2());
      CHECK(sc.schemaMajor == 2);
      CHECK(sc.keys.size() == 9);
      CHECK_NEAR(sc.cell90, 0.0838488197937161, 1e-12);
      const auto* dark = sc.Key(0., 180.);
      CHECK(dark != nullptr);
      if (dark) {
        CHECK(dark->dark1 && !dark->dark0);
        CHECK_NEAR(dark->tFlux0, 0.6, 0);
        CHECK_NEAR(dark->tFlux1, 0., 0);
        CHECK_NEAR(dark->cK, 1., 1e-12);
        CHECK_NEAR(dark->w1, 1., 0);
      }
      const auto* a = sc.Key(45., 135.);
      CHECK(a != nullptr && !a->dark0 && !a->dark1);
      if (a) { CHECK_REL(a->cK, std::cos(45. * CLHEP::pi / 180.), 1e-12); CHECK_NEAR(a->w0, 0.9, 0); CHECK_NEAR(a->tFlux1, 0.3, 0); }
      const auto* e = sc.Key(90., 90.);
      CHECK(e != nullptr && e->cK > 0.083 && e->cK < 0.085);
      CHECK(sc.Key(0., 100.) == nullptr && sc.Key(45. + 5e-10, 135.) != nullptr && sc.Key(45. + 2e-9, 135.) == nullptr);
      // declared order: phi-major, then theta
      CHECK(sc.keys[0].phiDeg == 0. && sc.keys[0].thetaDeg == 90. && sc.keys[1].thetaDeg == 135. && sc.keys[3].phiDeg == 45.);
      CHECK_NEAR(sc.exitY.step, 1e-3, 1e-15);
      CHECK_NEAR(sc.exitZ.step, 2.5e-5, 1e-18);
      CHECK(sc.exitFieldColumns.size() == 16 && sc.exitFieldColumns[15] == "InterceptFactor" && sc.exitFieldColumns[10] == "Ex_imag");
      CHECK(sc.lowestMode == "TE10" && !sc.invariant.empty());
      // The 1.x fixture is schema 1 with no per-key records.
      const auto s1 = Parse(Valid());
      CHECK(s1.schemaMajor == 1 && s1.keys.empty() && s1.cell90 == -1. && s1.exitY.step == 0.);
      // R2: a 2.x disc (uniform_disc, TE11 lowest) parses too.
      auto pd = hfssfix::NineKeyParams2(kStem);
      pd.base.disc = true; pd.base.radius = 5e-5;
      const auto disc = Parse(hfssfix::SidecarJson2(pd));
      CHECK(disc.schemaMajor == 2 && disc.lowestMode == "TE11" && disc.crossSection.shape == BBRCrossSection::kDisc);
      // Load reads a 2.x sidecar from disk as well.
      TempDir td;
      WriteFile(td.path() / "waveguides" / (kStem + ".dataset.json"), Valid2());
      CHECK(BBRDatasetSidecar::Load((td.path() / "waveguides").string(), "crack", kStem, 500.).keys.size() == 9);
    }},
    {"s2_version_gate", [] {
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["schema_version"] = "3.0"; })); }, "BBRDatasetSidecar",
                        "schema_version 3.0: BBRsim reads 1.x and 2.x");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["schema_version"] = "0.9"; })); }, "BBRDatasetSidecar",
                        "BBRsim reads 1.x and 2.x");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["schema_version"] = "20.0"; })); }, "BBRDatasetSidecar",
                        "BBRsim reads 1.x and 2.x");
      CHECK(Parse(Valid()).schemaMajor == 1);
      CHECK(Parse(Mut2([](json& j) { j["schema_version"] = "2.7"; j["x_future"] = {1, 2}; })).schemaMajor == 2);
      CHECK(Parse(Mut2([](json& j) { j["schema_version"] = "2"; })).schemaMajor == 2);
    }},
    {"s2_fixed_strings", [] {
      struct Edit1 { const char* pointer; json value; const char* field; };
      const std::vector<Edit1> edits = {
        {"/model/kind", "slot_in_wall", "model.kind"},
        {"/transmittance/definition", "OutgoingPower / IngoingPower", "transmittance.definition"},
        {"/transmittance/numerator", "exit_abs_re_s", "transmittance.numerator"},
        {"/transmittance/incoming", "area_times_intensity_times_cos", "transmittance.incoming"},
        {"/transmittance/above_one", "load_time_rescale", "transmittance.above_one"},
        {"/exit_field/tangential_components", json::array({"Ex", "Ez"}), "exit_field.tangential_components"},
        {"/exit_field/quadrature", "uniform_disc", "exit_field.quadrature"},
        {"/exit_field/sampler", "abs_e_squared", "exit_field.sampler"},
        {"/exit_field/scale/rule", "none", "exit_field.scale.rule"},
        {"/cross_term/density", "all_components", "cross_term.density"},
        {"/cross_term/above_limit", "ignore", "cross_term.above_limit"},
        {"/far_field/normalization", "absolute", "far_field.normalization"},
        {"/exit_field/column_order_check/rule", "none", "exit_field.column_order_check.rule"},
        {"/exit_field/column_order_check/reference/passed", false, "exit_field.column_order_check.reference.passed"},
        {"/exit_field/plane_tolerance_m", 1e-9, "exit_field.plane_tolerance_m"},
        {"/exit_field/columns", json::array({"Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower",
                                             "X", "Y", "Z", "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag",
                                             "Ez_imag", "InterceptFactor"}), "exit_field.columns"},
      };
      for (const auto& e : edits) {
        ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { At(j, e.pointer) = e.value; })); }, "BBRDatasetSidecar", e.field);
      }
      // trapezoid is the only quadrature of a rectangle, and uniform_disc of a disc.
      auto pd = hfssfix::NineKeyParams2(kStem);
      pd.base.disc = true; pd.base.radius = 5e-5;
      ExpectG4Exception("BBR025", [&] {
        json j = json::parse(hfssfix::SidecarJson2(pd));
        j["exit_field"]["quadrature"] = "trapezoid";
        Parse(j.dump());
      }, "BBRDatasetSidecar", "exit_field.quadrature");
      // A missing block is BBR024, a mistyped string too.
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j.erase("model"); })); }, "BBRDatasetSidecar", "missing field model");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["transmittance"]["numerator"] = 3; })); }, "BBRDatasetSidecar",
                        "transmittance.numerator must be a string");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"].erase("scale"); })); }, "BBRDatasetSidecar",
                        "missing field exit_field.scale");
    }},
    {"s2_per_key_lists", [] {
      struct List { const char* pointer; const char* number; bool pols; };   // number: a finite number inside element 0
      const std::vector<List> lists = {
        {"/transmittance/c_key", "c_k", false},
        {"/transmittance/per_key", "polarizations/0/t_flux", true},
        {"/exit_field/scale/per_key", "w0", false},
        {"/far_field/ff_over_flux", "polarizations/0", true},
      };
      for (const auto& l : lists) {
        const std::string what = l.pointer;
        const auto element = [&](json& j, std::size_t i) -> json& { return At(j, l.pointer)[i]; };
        // a duplicate key, a missing key, a key off the declared grid by 1e-6, a string angle
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { At(j, l.pointer).push_back(element(j, 0)); })); },
                          "BBRDatasetSidecar", "duplicate");
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { At(j, l.pointer).erase(At(j, l.pointer).size() - 1); })); },
                          "BBRDatasetSidecar", "no record for declared key");
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { element(j, 0)["theta_deg"] = 90. + 1e-6; })); },
                          "BBRDatasetSidecar", "does not match a declared key");
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { element(j, 0)["phi_deg"] = "0"; })); },
                          "BBRDatasetSidecar", "phi_deg must be a finite number");
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { At(j, l.pointer) = json::object(); })); },
                          "BBRDatasetSidecar", "must be a list");
        if (l.pols) {
          ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { element(j, 0)["polarizations"].erase("1"); })); },
                            "BBRDatasetSidecar", "polarizations");
          ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { element(j, 0).erase("polarizations"); })); },
                            "BBRDatasetSidecar", "polarizations");
        }
        // a non-finite number: 1e999 and NaN stop in the JSON parser, a string stands in the field check
        ExpectG4Exception("BBR024", [&] {
          json j = json::parse(Valid2());
          element(j, 0)[json::json_pointer("/" + std::string(l.number))] = 12345.6789;
          Parse(Sentinel(j.dump(2), "1e999"));
        }, "BBRDatasetSidecar", "number overflow");
        ExpectG4Exception("BBR024", [&] {
          json j = json::parse(Valid2());
          element(j, 0)[json::json_pointer("/" + std::string(l.number))] = 12345.6789;
          Parse(Sentinel(j.dump(2), "NaN"));
        }, "BBRDatasetSidecar", "not valid JSON");
        ExpectG4Exception("BBR024", [&] { Parse(Mut2([&](json& j) { element(j, 0)[json::json_pointer("/" + std::string(l.number))] = "nan"; })); },
                          "BBRDatasetSidecar", "finite number");
        // 5e-10 off the declared key still matches (1e-9)
        Parse(Mut2([&](json& j) { element(j, 0)["theta_deg"] = 90. + 5e-10; element(j, 3)["phi_deg"] = 45. - 5e-10; }));
      }
      // ff_over_flux: null on a lit basis is refused, null on a dark basis parses (the fixture has one)
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["far_field"]["ff_over_flux"][0]["polarizations"]["0"] = nullptr; })); },
                        "BBRDatasetSidecar", "null");
      CHECK(json::parse(Valid2())["far_field"]["ff_over_flux"][2]["polarizations"]["1"].is_null());
      Parse(Valid2());
      // a scale weight must be positive
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["scale"]["per_key"][0]["w0"] = 0.; })); },
                        "BBRDatasetSidecar", "w0");
      // t_flux and c_k of per_key must agree with c_key
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["transmittance"]["per_key"][1]["c_k"] = 0.5; })); },
                        "BBRDatasetSidecar", "c_k");
    }},
    {"s2_c_k_and_cell90", [] {
      CHECK_NEAR(BBRDatasetSidecar::Cell90({0, 45, 90, 135, 180}), 0.08384881979, 1e-10);
      CHECK_NEAR(BBRDatasetSidecar::Cell90({90, 105, 120, 135, 150, 165, 180}), 0.03663305469, 1e-10);
      CHECK_NEAR(BBRDatasetSidecar::Cell90({0, 45, 90, 135, 180}), 0.0838488197937161, 1e-14);
      CHECK_NEAR(BBRDatasetSidecar::Cell90({90, 105, 120, 135, 150, 165, 180}), 0.036633054690581604, 1e-14);
      CHECK(BBRDatasetSidecar::Cell90({120, 180}) == -1.);
      CHECK(BBRDatasetSidecar::Cell90({0, 45, 90}) == -1.);   // nothing above 90 deg
      CHECK_REL(BBRDatasetSidecar::InterceptFactor(135., {90, 135, 180}), std::cos(45. * CLHEP::pi / 180.), 1e-14);
      CHECK_NEAR(BBRDatasetSidecar::InterceptFactor(180., {90, 135, 180}), 1., 1e-15);
      CHECK_NEAR(BBRDatasetSidecar::InterceptFactor(90., {90, 135, 180}), BBRDatasetSidecar::Cell90({90, 135, 180}), 0);
      CHECK_REL(BBRDatasetSidecar::InterceptFactor(45., {0, 45, 90, 135, 180}), std::cos(45. * CLHEP::pi / 180.), 1e-14);
      // a declared grid without 90 deg cannot give the cell: refused
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["excitation"]["incident_theta_deg"] = {135., 180.}; })); },
                        "BBRDatasetSidecar", "theta = 90");
      // c_key at (0, 135) off by 2e-6 relative is refused (BBR025), off by 5e-7 parses
      const double c135 = std::cos(45. * CLHEP::pi / 180.);
      ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { j["transmittance"]["c_key"][1]["c_k"] = c135 * (1. + 2e-6); })); },
                        "BBRDatasetSidecar", "c_k");
      Parse(Mut2([&](json& j) { j["transmittance"]["c_key"][1]["c_k"] = c135 * (1. + 5e-7); }));
      // the cell: cell_90 and the c_k at 90 deg must follow the declared grid
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90"] = 0.08384881979 * 1.001; })); },
                        "BBRDatasetSidecar", "cell_90");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["transmittance"]["c_key"][0]["c_k"] = 0.19134; })); },
                        "BBRDatasetSidecar", "c_k");
      Parse(Mut2([](json& j) { j["transmittance"]["cell_90"] = 0.0838488197937161 * (1. + 5e-7); }));
      // the 35-key style grid: Delta = 15 deg gives the 0.0366 cell
      auto p = hfssfix::NineKeyParams2(kStem);
      p.base.thetas = {90., 105., 120., 135., 150., 165., 180.};
      p.tflux.clear();
      const auto sc = Parse(hfssfix::SidecarJson2(p));
      CHECK_NEAR(sc.cell90, 0.03663305469, 1e-10);
      CHECK(sc.keys.size() == 21);
    }},
    {"s2_dark_tflux", [] {
      // OR-038: a dark basis has t_flux exactly 0; the sidecar contradicting itself is BBR024.
      const auto bad = Mut2([](json& j) { j["transmittance"]["per_key"][2]["polarizations"]["1"]["t_flux"] = 1e-15; });
      ExpectG4Exception("BBR024", [&] { Parse(bad); }, "BBRDatasetSidecar", "key (phi 0, theta 180)");
      ExpectG4Exception("BBR024", [&] { Parse(bad); }, "BBRDatasetSidecar", "polarization 1");
      ExpectG4Exception("BBR024", [&] { Parse(bad); }, "BBRDatasetSidecar", "1e-15");
      ExpectG4Exception("BBR024", [&] { Parse(bad); }, "BBRDatasetSidecar", "OR-038");
      // dark0 on a lit key with its t_flux
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["scale"]["per_key"][0]["dark0"] = true; })); },
                        "BBRDatasetSidecar", "polarization 0");
      // dark and lit with t_flux 0 is fine; the fixture's own dark key parses
      Parse(Valid2());
    }},
    {"s2_emitter", [] {
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"]["emitter"] = "lambertian"; })); },
                        "BBRDatasetSidecar", "cell_90_assumes.emitter");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"]["cos_floor"] = 0.02; })); },
                        "BBRDatasetSidecar", "cell_90_assumes.cos_floor");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"]["band_deg"] = {60., 90.}; })); },
                        "BBRDatasetSidecar", "cell_90_assumes.band_deg");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"].erase("band_deg"); })); },
                        "BBRDatasetSidecar", "band_deg");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"]["band_deg"] = {67.5}; })); },
                        "BBRDatasetSidecar", "band_deg");
      Parse(Mut2([](json& j) { j["transmittance"]["cell_90_assumes"]["band_deg"] = {67.5 + 5e-10, 90.}; }));
    }},
    {"s2_grid_step", [] {
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["y_e"]["step"] = 1.001e-3; })); },
                        "BBRDatasetSidecar", "exit_field.grid.y_e.step");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["z_e"]["step"] = 2.4e-5; })); },
                        "BBRDatasetSidecar", "exit_field.grid.z_e.step");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["y_e"].erase("step"); })); },
                        "BBRDatasetSidecar", "missing field exit_field.grid.y_e.step");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["z_e"].erase("step"); })); },
                        "BBRDatasetSidecar", "exit_field.grid.z_e.step");
      // Choice 3: a rectangle axis with one node has no step and no trapezoid
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["y_e"]["count"] = 1; })); },
                        "BBRDatasetSidecar", "exit_field.grid.y_e needs count >= 2");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["exit_field"]["grid"]["z_e"]["count"] = 1; })); },
                        "BBRDatasetSidecar", "exit_field.grid.z_e needs count >= 2");
      // 1e-6 relative on the step is within tolerance
      Parse(Mut2([](json& j) { j["exit_field"]["grid"]["y_e"]["step"] = 1e-3 * (1. + 5e-7); }));
      // The 1.x grid has no step and does not need one.
      CHECK(Parse(Valid()).exitY.step == 0.);
    }},
    {"s2_regime_cutoffs", [] {
      using RC = std::vector<std::pair<std::string, double>>;
      BBRCrossSection slot1;   // 10 mm x 50 um: crack1
      slot1.shape = BBRCrossSection::kRectangle; slot1.yHalf_m = 5e-3; slot1.zHalf_m = 25e-6;
      const RC r1 = BBRDatasetSidecar::RegimeCutoffs(slot1);
      CHECK(r1.size() == 7);
      if (r1.size() == 7) {
        CHECK(r1[0].first == "TE10" && r1[1].first == "TE01" && r1[6].first == "TE06");
        CHECK_REL(r1[0].second, 14.9896229, 1e-6);
        CHECK_REL(r1[1].second, 2997.92458, 1e-6);
        CHECK_REL(r1[6].second, 17987.5475, 1e-6);
      }
      BBRCrossSection slot2 = slot1;   // 10 mm x 100 um: crack2
      slot2.zHalf_m = 50e-6;
      const RC r2 = BBRDatasetSidecar::RegimeCutoffs(slot2);
      CHECK(r2.size() == 14);
      if (r2.size() == 14) {
        CHECK(r2[0].first == "TE10" && r2[9].first == "TE09" && r2[10].first == "TE0,10" && r2[13].first == "TE0,13");
        CHECK_REL(r2[13].second, 19486.5098, 1e-6);
      }
      BBRCrossSection disc;
      disc.shape = BBRCrossSection::kDisc; disc.radius_m = 50e-6;
      const RC rd = BBRDatasetSidecar::RegimeCutoffs(disc);
      CHECK(rd.size() == 3);
      if (rd.size() == 3) {
        CHECK(rd[0].first == "TE11" && rd[1].first == "TM01" && rd[2].first == "TE21");
        CHECK_REL(rd[0].second, 1756.98, 1e-5);
        CHECK_REL(rd[1].second, 2294.85, 1e-5);
        CHECK_REL(rd[2].second, 2914.56, 1e-5);
      }
      // The fixture's own list (rectangle and disc) is accepted and parsed; 1.x gives the lowest mode only.
      const auto sc = Parse(Valid2());
      CHECK(sc.regimeCutoffs.size() == 7 && sc.regimeCutoffs[0].first == "TE10");
      CHECK_REL(sc.regimeCutoffs[0].second, sc.cutoffGHz, 1e-12);
      auto pd = hfssfix::NineKeyParams2(kStem);
      pd.base.disc = true; pd.base.radius = 5e-5;
      CHECK(Parse(hfssfix::SidecarJson2(pd)).regimeCutoffs.size() == 3);
      const auto s1 = Parse(Valid());
      CHECK(s1.regimeCutoffs.size() == 1 && s1.regimeCutoffs[0].first == s1.lowestMode &&
            s1.regimeCutoffs[0].second == s1.cutoffGHz && s1.rhoValidBelowGHz == -1.);
      // An entry off by 2e-6 relative is refused (BBR025), one off by 5e-7 parses.
      const double c1 = sc.regimeCutoffs[1].second;
      ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { j["modes"]["regime_cutoffs_ghz"][1]["cutoff_ghz"] = c1 * (1. + 2e-6); })); },
                        "BBRDatasetSidecar", "regime_cutoffs_ghz");
      ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { j["modes"]["regime_cutoffs_ghz"][1]["cutoff_ghz"] = c1 * (1. + 2e-6); })); },
                        "BBRDatasetSidecar", "TE01");
      Parse(Mut2([&](json& j) { j["modes"]["regime_cutoffs_ghz"][1]["cutoff_ghz"] = c1 * (1. + 5e-7); }));
      // A missing, an extra and a relabelled entry are refused.
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["modes"]["regime_cutoffs_ghz"].erase(6); })); },
                        "BBRDatasetSidecar", "regime_cutoffs_ghz");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["modes"]["regime_cutoffs_ghz"].push_back({{"mode", "TE07"}, {"cutoff_ghz", 20983.}}); })); },
                        "BBRDatasetSidecar", "regime_cutoffs_ghz");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& j) { j["modes"]["regime_cutoffs_ghz"][2]["mode"] = "TE0,2"; })); },
                        "BBRDatasetSidecar", "TE0,2");
      // Mistyped or absent: BBR024.
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["modes"].erase("regime_cutoffs_ghz"); })); },
                        "BBRDatasetSidecar", "missing field modes.regime_cutoffs_ghz");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["modes"]["regime_cutoffs_ghz"] = "TE10"; })); },
                        "BBRDatasetSidecar", "modes.regime_cutoffs_ghz must be an array");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["modes"]["regime_cutoffs_ghz"][0].erase("mode"); })); },
                        "BBRDatasetSidecar", "mode");
      // A section taller along g than l: TE01 is the lowest, so the ascending list starts with it.
      auto pt = hfssfix::NineKeyParams2(kStem);
      pt.base.yHalf = 1e-5; pt.base.zHalf = 2e-5;
      const auto tall = Parse(hfssfix::SidecarJson2(pt));
      CHECK(!tall.regimeCutoffs.empty() && tall.regimeCutoffs[0].first == "TE01");
      for (std::size_t i = 1; i < tall.regimeCutoffs.size(); ++i) CHECK(tall.regimeCutoffs[i].second >= tall.regimeCutoffs[i - 1].second);
    }},
    {"s2_dropped_modes_fields", [] {
      // Schema 2 drops the lists, counts, list_limit_ghz and basis (OR-031): a block without them parses.
      const auto plain = Parse(Mut2([](json& j) {
        for (const char* k : {"propagating_count", "basis", "list_limit_ghz", "gap_family_onsets", "mode_count_below_limit", "cutoffs"})
          j["modes"].erase(k);
      }));
      CHECK(plain.schemaMajor == 2 && plain.regimeCutoffs.size() == 7);
      // One that still carries the 1.x lists parses too: unknown fields are ignored, never compared.
      const auto legacy = Parse(Mut2([](json& j) {
        j["modes"]["propagating_count"] = 33;
        j["modes"]["basis"] = "closed PEC a=0.009 b=5e-05";
        j["modes"]["list_limit_ghz"] = 20000.;
        j["modes"]["gap_family_onsets"] = json::array({{{"n", 1}, {"mode", "TE01"}, {"cutoff_ghz", 1.}}});
        j["modes"]["mode_count_below_limit"] = 9999;
        j["modes"]["cutoffs"] = json::array({{{"mode", "TE11"}, {"cutoff_ghz", 1.}, {"degeneracy", 2}}});
      }));
      CHECK(legacy.regimeCutoffs.size() == 7 && legacy.SameInvariant(plain));
    }},
    {"s2_cross_term", [] {
      const auto sc = Parse(Valid2());
      CHECK_REL(sc.rhoValidBelowGHz, sc.polarizationFilterLimitGHz, 1e-12);
      CHECK_REL(sc.rhoValidBelowGHz, 2997.92458, 1e-6);   // c/(2b), b = 50 um
      const double lim = sc.polarizationFilterLimitGHz;
      ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { j["cross_term"]["valid_below_ghz"] = lim * (1. + 2e-6); })); },
                        "BBRDatasetSidecar", "cross_term.valid_below_ghz");
      ExpectG4Exception("BBR025", [&] { Parse(Mut2([&](json& j) { j["cross_term"]["valid_below_ghz"] = lim * (1. - 2e-6); })); },
                        "BBRDatasetSidecar", "polarization_filter_limit_ghz");
      const auto near = Parse(Mut2([&](json& j) { j["cross_term"]["valid_below_ghz"] = lim * (1. + 5e-7); }));
      CHECK_REL(near.rhoValidBelowGHz, lim * (1. + 5e-7), 1e-12);
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& j) { j["cross_term"].erase("valid_below_ghz"); })); },
                        "BBRDatasetSidecar", "cross_term.valid_below_ghz");
    }},
    {"invariant_schema_and_grid", [] {
      using V = std::vector<std::string>;
      // Choice 1: one ID must not mix schema majors or incidence grids across frequencies.
      const auto one = Parse(Valid());
      const auto two = Parse(Valid2());
      const V d = one.InvariantDiff(two);
      CHECK(!d.empty() && d[0] == "schema" && !one.SameInvariant(two) && !two.SameInvariant(one));
      CHECK(two.InvariantDiff(one) == one.InvariantDiff(two));
      // A 15-key style and a 35-key style grid of one geometry differ in "excitation" only.
      auto p15 = hfssfix::NineKeyParams2(kStem);
      auto p35 = hfssfix::NineKeyParams2(kStem);
      p35.base.thetas = {90., 105., 120., 135., 150., 165., 180.};
      p35.base.phis = {0., 22.5, 45., 67.5, 90.};
      p35.tflux.clear();
      const auto a = Parse(hfssfix::SidecarJson2(p15)), b = Parse(hfssfix::SidecarJson2(p35));
      CHECK(a.InvariantDiff(b) == V{"excitation"} && b.InvariantDiff(a) == V{"excitation"} && !a.SameInvariant(b));
      // The same grid at another frequency agrees (2.x), and so do two 1.x sidecars of one grid (R3).
      const auto a2 = Parse(hfssfix::SidecarJson2(p15));
      CHECK(a.InvariantDiff(a2).empty() && a.SameInvariant(a2));
      const std::string stem300 = "crack_300GHz";
      const auto s300 = Parse(hfssfix::SidecarJson(hfssfix::SidecarFrom(stem300, hfssfix::Mini500())), stem300, 300.);
      CHECK(one.InvariantDiff(s300).empty() && one.SameInvariant(s300) && s300.SameInvariant(one));
      // A changed phi list alone is also a different grid.
      auto pphi = hfssfix::NineKeyParams2(kStem);
      pphi.base.phis = {0., 90.};
      pphi.tflux.clear();
      CHECK(a.InvariantDiff(Parse(hfssfix::SidecarJson2(pphi))) == V{"excitation"});
      // Names come in order: schema first, then the existing blocks.
      auto pm = hfssfix::NineKeyParams2(kStem);
      pm.base.yHalf = 4.0e-3;
      const auto sec = Parse(hfssfix::SidecarJson2(pm));
      CHECK(one.InvariantDiff(sec).front() == "schema");
    }},
    {"f12_fits_solid", [] {
      const auto box = Parse(Valid());                  // section 4.5 mm x 25 um
      box.CheckFitsSolid(G4Box("ok", 2 * mm, 4.5 * mm, 0.025 * mm), "ok");
      // The fit message names the volume; the origin message is the other check of CheckFitsSolid.
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("thin", 2 * mm, 4.5 * mm, 0.024 * mm), "thin"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not equal crack volume thin");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("short", 2 * mm, 4 * mm, 0.025 * mm), "short"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not equal crack volume short");
      // A solid whose local origin is not inside it: the wrapper finds the exit face from the origin.
      G4DisplacedSolid off("off", new G4Box("b", 2 * mm, 5 * mm, 0.026 * mm), G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 1 * mm)));
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(off, "off"); }, "BBRDatasetSidecar",
                        "the local origin of crack volume off is not inside its solid");
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5; p.lengthMm = 0.4;   // Tube() is 0.4 mm long
      const auto disc = Parse(hfssfix::SidecarJson(p));
      disc.CheckFitsSolid(*Tube(0.050 * mm), "tube50");
      // A plain G4Tubs keeps its axis on local z, so local x is radial: rejected at startup.
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(G4Tubs("bare", 0., 0.050 * mm, 0.2 * mm, 0., CLHEP::twopi), "bare"); }, "BBRDatasetSidecar",
                        "the declared exit cross-section does not equal crack volume bare");
    }},
    {"f12_length", [] {
      // CDX-028, OR-020: the crack volume's length (DistanceToOut along +x and -x from the origin) must
      // equal geometry.extent_mm.p within 1 nm, or 1 nm plus 1e-6 relative where longer.
      const auto box = Parse(Valid());                  // p = 4 mm, section 9 mm x 50 um
      box.CheckFitsSolid(G4Box("ok", 2 * mm, 4.5 * mm, 0.025 * mm), "ok");
      box.CheckFitsSolid(G4Box("longer_half_nm", 2 * mm + 5e-7 * mm, 4.5 * mm, 0.025 * mm), "longer_half_nm");
      box.CheckFitsSolid(G4Box("shorter_half_nm", 2 * mm - 4e-7 * mm, 4.5 * mm, 0.025 * mm), "shorter_half_nm");
      // Longer: 1 nm plus 1e-6 relative of the 4 mm, so L + 4e-6 mm (half-lengths + 2e-6 mm) is inside, L + 7e-6 mm is not.
      box.CheckFitsSolid(G4Box("longer_4nm_total", 2 * mm + 2e-6 * mm, 4.5 * mm, 0.025 * mm), "longer_4nm_total");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("longer_7nm_total", 2 * mm + 3.5e-6 * mm, 4.5 * mm, 0.025 * mm), "longer_7nm_total"); },
                        "BBRDatasetSidecar", "length");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("shorter_2nm", 2 * mm - 2e-6 * mm, 4.5 * mm, 0.025 * mm), "shorter_2nm"); },
                        "BBRDatasetSidecar", "length");
      // The message names both lengths and OR-020.
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("slab", 2.5 * mm, 4.5 * mm, 0.025 * mm), "slab"); },
                        "BBRDatasetSidecar", "crack volume's length 5 mm differs from the HFSS model's 4 mm (OR-020)");
      // A disc: a 0.4 mm declared length against the 0.4 mm tube passes, a 0.5 mm tube does not.
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5; p.lengthMm = 0.4;
      const auto disc = Parse(hfssfix::SidecarJson(p));
      disc.CheckFitsSolid(*Tube(0.050 * mm), "tube50");
      G4RotationMatrix ry;
      ry.rotateY(90. * deg);
      G4DisplacedSolid longTube("long_tube", new G4Tubs("t", 0., 0.050 * mm, 0.25 * mm, 0., CLHEP::twopi),
                                G4Transform3D(ry, G4ThreeVector()));
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(longTube, "long_tube"); }, "BBRDatasetSidecar", "length");
    }},
    {"f11_extent_section", [] {
      // CDX-028: a rectangle's extent_mm.l and g must equal the full widths of its own declared section.
      auto with = [](const std::string& key, double v) {
        json j = json::parse(Valid());
        j["geometry"]["extent_mm"][key] = v;
        return j.dump(2);
      };
      Parse(Valid());                                   // l = 9 mm, g = 0.05 mm: consistent
      ExpectG4Exception("BBR024", [&] { Parse(with("l", 9.00001)); }, "BBRDatasetSidecar", "geometry.extent_mm.l");
      ExpectG4Exception("BBR024", [&] { Parse(with("l", 9.00001)); }, "BBRDatasetSidecar", "exit_field.cross_section");
      ExpectG4Exception("BBR024", [&] { Parse(with("g", 0.05 * (1 + 2e-6))); }, "BBRDatasetSidecar", "geometry.extent_mm.g");
      Parse(with("g", 0.05 * (1 + 5e-7)));              // inside 1e-6 relative
      Parse(with("l", 9. * (1 - 5e-7)));
      // A disc's extent is not compared with its radius here; the round gap's extent_mm is the HFSS bounding box.
      json j = json::parse(Valid());
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5; p.lengthMm = 0.4;
      json d = json::parse(hfssfix::SidecarJson(p));
      d["geometry"]["extent_mm"]["l"] = 0.1234;
      Parse(d.dump(2));
    }},
    {"material_vacuum", [] {
      // CDX-029: provenance.material, when a string, must be "vacuum" (case-insensitive).
      auto with = [](const json& v) {
        json j = json::parse(Valid());
        j["provenance"]["material"] = v;
        return j.dump(2);
      };
      Parse(with("vacuum"));
      Parse(with("Vacuum"));
      Parse(with("VACUUM"));
      ExpectG4Exception("BBR025", [&] { Parse(with("silicon")); }, "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR025", [&] { Parse(with("silicon")); }, "BBRDatasetSidecar", "silicon");
      ExpectG4Exception("BBR025", [&] { Parse(with("")); }, "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR024", [&] { Parse(with(3)); }, "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR024", [&] { Parse(with(true)); }, "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR024", [&] { Parse(with(json::array({"vacuum"}))); }, "BBRDatasetSidecar", "provenance.material");
    }},
    {"material_null_schema1", [] {
      // The legacy 500 GHz sidecars carry "material": null; null and an absent key stay accepted under 1.x.
      json j = json::parse(Valid());
      j["provenance"]["material"] = nullptr;
      Parse(j.dump(2));
      j["provenance"].erase("material");
      Parse(j.dump(2));
      j.erase("provenance");                            // no provenance object at all
      Parse(j.dump(2));
    }},
    {"material_required_schema2", [] {
      // CDX-029, OR-039: under 2.x provenance.material is required; null, an absent key or an absent
      // provenance is BBR024 naming the field. A string still follows the shared rule.
      json j = json::parse(Valid2());
      CHECK(j.at("provenance").at("material") == "vacuum");   // the 2.0 fixture writes it
      j["provenance"]["material"] = "Vacuum";
      Parse(j.dump(2));
      j["provenance"]["material"] = nullptr;
      ExpectG4Exception("BBR024", [&] { Parse(j.dump(2)); }, "BBRDatasetSidecar", "provenance.material");
      j["provenance"].erase("material");
      ExpectG4Exception("BBR024", [&] { Parse(j.dump(2)); }, "BBRDatasetSidecar", "provenance.material");
      j.erase("provenance");
      ExpectG4Exception("BBR024", [&] { Parse(j.dump(2)); }, "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR025", [] { Parse(Mut2([](json& m) { m["provenance"]["material"] = "silicon"; })); },
                        "BBRDatasetSidecar", "provenance.material");
      ExpectG4Exception("BBR024", [] { Parse(Mut2([](json& m) { m["provenance"]["material"] = 3; })); },
                        "BBRDatasetSidecar", "provenance.material");
    }},
    {"f12_equality", [] {
      // OR-020: the solid must equal the declared section, so a larger one fails as well as a smaller one.
      const auto box = Parse(Valid());                  // section 4.5 mm x 25 um
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("wide", 2 * mm, 4.5 * mm, 0.026 * mm), "wide"); }, "BBRDatasetSidecar",
                        "larger than");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("narrow", 2 * mm, 4.5 * mm, 0.024 * mm), "narrow"); }, "BBRDatasetSidecar",
                        "smaller than");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("long", 2 * mm, 5 * mm, 0.025 * mm), "long"); }, "BBRDatasetSidecar",
                        "larger than");
      // Tolerance edges on the 25 um half-gap: 1 nm absolute smaller, 1 nm (plus 1e-6 relative) larger.
      box.CheckFitsSolid(G4Box("larger_half_nm", 2 * mm, 4.5 * mm, 0.025 * mm + 5e-7 * mm), "larger_half_nm");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("larger_2nm", 2 * mm, 4.5 * mm, 0.025 * mm + 2e-6 * mm), "larger_2nm"); },
                        "BBRDatasetSidecar", "larger than");
      ExpectG4Exception("BBR025", [&] { box.CheckFitsSolid(G4Box("smaller_2nm", 2 * mm, 4.5 * mm, 0.025 * mm - 2e-6 * mm), "smaller_2nm"); },
                        "BBRDatasetSidecar", "smaller than");
      auto p = hfssfix::SidecarFrom(kStem, hfssfix::Mini500());
      p.disc = true; p.radius = 5e-5; p.lengthMm = 0.4;
      const auto disc = Parse(hfssfix::SidecarJson(p));
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(*Tube(0.051 * mm), "tube51"); }, "BBRDatasetSidecar", "larger than");
      ExpectG4Exception("BBR025", [&] { disc.CheckFitsSolid(*Tube(0.049 * mm), "tube49"); }, "BBRDatasetSidecar", "smaller than");
    }},
  });
}
