#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "BBRConfigManager.hh"
#include "BBRCrackLibrary.hh"
#include "BBRMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <nlohmann/json.hpp>

#include <cmath>
#include <sstream>
#include <set>
#include <thread>
#include <vector>

using namespace bbrtest;

namespace {
// A grid of frequencies for one id: writes Mini500 under <id>_<tok>GHz for each token.
void Grid(const std::filesystem::path& root, const std::string& id, std::initializer_list<const char*> toks) {
  for (const char* t : toks) {
    auto ds = hfssfix::Mini500();
    // the Freq column must agree with the directory (BBR009): rewrite it
    const std::string f = std::string(t) + "GHz";
    // resume the search after each replacement: f may itself contain "500GHz" (500, 1500)
    auto fix = [&](std::string s) {
      for (std::size_t p = 0; (p = s.find("500GHz", p)) != std::string::npos; p += f.size()) s.replace(p, 6, f);
      return s;
    };
    ds.ff0 = fix(ds.ff0); ds.ff1 = fix(ds.ff1); ds.wg0 = fix(ds.wg0); ds.wg1 = fix(ds.wg1);
    hfssfix::WriteDataset(root, id + "_" + f, ds);
  }
}
// A dataset on NineKeyGrid2 (schema 2) or NineKeyGrid (schema 1) at one frequency: a crack2-like
// rectangle (a = 10 mm, b = 100 um) or, for `disc`, a 50 um disc whose exit lattice is scaled into it.
void GridModes(const std::filesystem::path& root, const std::string& id, int freq, bool disc, bool schema1 = false) {
  auto d = schema1 ? hfssfix::NineKeyGrid() : hfssfix::NineKeyGrid2();
  const std::string f = std::to_string(freq) + "GHz";
  auto fix = [&](const std::string& csv, bool scaleY) {
    std::istringstream in(csv);
    std::string line, out;
    bool header = true;
    while (std::getline(in, line)) {
      if (!header) {
        for (std::size_t p = 0; (p = line.find("500GHz,", p)) != std::string::npos; p += f.size()) line.replace(p, 6, f);
        if (scaleY) {
          std::vector<std::string> v;
          std::stringstream ss(line);
          std::string t;
          while (std::getline(ss, t, ',')) v.push_back(t);
          if (v.size() > 8) {   // a waveguide row; far-field rows have 10 fields and are left alone
            v[7] = hfssfix::detail::Fmt(std::strtod(v[7].c_str(), nullptr) * 2e-2);
            line.clear();
            for (std::size_t i = 0; i < v.size(); ++i) line += (i ? "," : "") + v[i];
          }
        }
      }
      header = false;
      out += line + "\n";
    }
    return out;
  };
  d.ff0 = fix(d.ff0, false); d.ff1 = fix(d.ff1, false);
  d.wg0 = fix(d.wg0, disc); d.wg1 = fix(d.wg1, disc);
  const std::string stem = id + "_" + f;
  hfssfix::Schema2Params p = hfssfix::Schema2From(stem, d);
  p.tflux.assign(p.base.phis.size() * p.base.thetas.size(), std::array<double, 4>{0.6, 0.3, 0., 0.});
  p.tflux[2] = {0.6, 0., 0., 1.};
  if (disc) { p.base.disc = true; p.base.radius = 5e-5; }
  else { p.base.yHalf = 5e-3; p.base.zHalf = 5e-5; }
  if (schema1) {
    const auto wg = hfssfix::WriteDataset(root, stem, d, false);
    bbrtest::WriteFile(wg / (stem + ".dataset.json"), hfssfix::SidecarJson(p.base));
  } else {
    hfssfix::WriteDataset2(root, stem, d, p);
  }
}
double Pick(const char* id, double nu) { double chosen = -1; BBRCrackLibrary::Instance().Lookup(id, nu, chosen); return chosen; }
}  // namespace

int main(int argc, char** argv) {
  return RunCase(argc, argv, {
    {"selection_grid5", [] {
      TempDir d; Grid(d.path(), "g5", {"50", "150", "500", "1500", "5000"});
      Grid(d.path(), "one", {"500"}); Grid(d.path(), "lin", {"100", "1000"}); Grid(d.path(), "tie", {"1", "100"});
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK_NEAR(Pick("g5", 500), 500, 0);
      // +-1 % around each log midpoint picks the nearer side in log space
      for (auto [lo, hi] : {std::pair{50., 150.}, {150., 500.}, {500., 1500.}, {1500., 5000.}}) {
        const double mid = std::sqrt(lo * hi);
        CHECK_NEAR(Pick("g5", mid * 0.99), lo, 0);
        CHECK_NEAR(Pick("g5", mid * 1.01), hi, 0);
      }
      CHECK_NEAR(Pick("g5", 50), 50, 0);   CHECK_NEAR(Pick("g5", 5000), 5000, 0);   // exact edges, on grid
      CHECK(Handler().CountWarnings("BBR008") == 0);
      CHECK_NEAR(Pick("g5", 20), 50, 0);   CHECK_NEAR(Pick("g5", 10), 50, 0);       // clamp low
      CHECK_NEAR(Pick("g5", 1e4), 5000, 0); CHECK_NEAR(Pick("g5", 2e4), 5000, 0);   // clamp high
      CHECK(Handler().CountWarnings("BBR008") == 2);                                  // once per side
      CHECK_NEAR(Pick("tie", 10), 1, 0);            // exact log tie -> lower
      CHECK_NEAR(Pick("lin", 400), 1000, 0);        // log-nearest, not linear-nearest
      CHECK_NEAR(Pick("lin", 200), 100, 0);
      CHECK_NEAR(Pick("one", 1), 500, 0); CHECK_NEAR(Pick("one", 1e7), 500, 0);
      CHECK(Handler().CountWarnings("BBR008") == 2);   // a one-point grid never warns
      // The fixture's TE10 cutoff is 16.66 GHz. A photon clamped across it on a
      // multi-point grid raises both BBR008 and BBR026 (g5's low clamps to 50 GHz:
      // BBR008 at 20 GHz, BBR026 at 10 GHz, the one below the cutoff), and a
      // one-point grid raises BBR026 though it never clamps (one at 1 GHz, served
      // by 500 GHz). One per dataset and direction: 2.
      CHECK(Handler().CountWarnings("BBR026") == 2);
    }},
    {"token_verbatim", [] {
      TempDir d; Grid(d.path(), "sci", {"1.5e3"});
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK_NEAR(Pick("sci", 1000), 1500, 0);
    }},
    {"discovery_errors", [] {
      TempDir d;
      hfssfix::WriteDataset(d.path(), "legacy", hfssfix::Mini500());                 // legacy <id>_Ephi=0 (no frequency)
      Grid(d.path(), "legacy", {"500"});             // + a valid dir beside it: only the legacy-specific BBR011 trigger explains a fatal here
      Grid(d.path(), "dup", {"500", "5e2"});
      std::filesystem::create_directories(d.path() / "waveguides" / "junk_abcGHz_Ephi=0");
      Grid(d.path(), "half", {"500"});
      std::filesystem::remove_all(d.path() / "waveguides" / "half_500GHz_Ephi=1");
      BBRConfigManager::SetDataDir(d.path().string());
      ExpectG4Exception("BBR011", [] { Pick("legacy", 500); }, "BBRCrackLibrary");
      ExpectG4Exception("BBR011", [] { Pick("dup", 500); }, "BBRCrackLibrary");
      ExpectG4Exception("BBR011", [] { Pick("junk", 500); }, "BBRCrackLibrary");
      ExpectG4Exception("BBR011", [] { Pick("nosuch", 500); }, "BBRCrackLibrary");
      // Ephi=1 missing: at discovery, before any CSV is read, and again (nothing was cached, the lock was released)
      ExpectG4Exception("BBR001", [] { Pick("half", 500); }, "BBRCrackLibrary::Discover", "half_500GHz_Ephi=1/far_field.csv");
      ExpectG4Exception("BBR001", [] { Pick("half", 500); }, "BBRCrackLibrary::Discover", "half_500GHz_Ephi=1/far_field.csv");
    }},
    {"discovery_missing_far_field", [] {
      // One CSV of four missing, at the frequency a photon would not select: BBR001 at discovery.
      TempDir d; Grid(d.path(), "ff", {"50", "500"});
      const auto gone = d.path() / "waveguides" / "ff_500GHz_Ephi=0" / "far_field.csv";
      std::filesystem::remove(gone);
      BBRConfigManager::SetDataDir(d.path().string());
      ExpectG4Exception("BBR001", [] { Pick("ff", 50); }, "BBRCrackLibrary::Discover",
                        "Missing HFSS file " + std::filesystem::absolute(gone).string() +
                          " (no such regular file): every frequency of dataset ff needs far_field.csv and "
                          "waveguide.csv in both ff_500GHz_Ephi=0 and ff_500GHz_Ephi=1.");
    }},
    {"discovery_missing_waveguide", [] {
      // A waveguide.csv that is a directory is no regular file either: BBR002 at discovery.
      TempDir d; Grid(d.path(), "wg", {"500"});
      const auto bad = d.path() / "waveguides" / "wg_500GHz_Ephi=1" / "waveguide.csv";
      std::filesystem::remove(bad);
      std::filesystem::create_directories(bad);
      BBRConfigManager::SetDataDir(d.path().string());
      ExpectG4Exception("BBR002", [] { Pick("wg", 500); }, "BBRCrackLibrary::Discover",
                        "Missing HFSS file " + std::filesystem::absolute(bad).string() +
                          " (no such regular file): every frequency of dataset wg needs far_field.csv and "
                          "waveguide.csv in both wg_500GHz_Ephi=0 and wg_500GHz_Ephi=1.");
    }},
    {"dataset_id_of", [] {
      CHECK(BBRCrackLibrary::DatasetIdOf("gap:1") == "gap");
      CHECK(BBRCrackLibrary::DatasetIdOf("gap") == "gap");
      CHECK(BBRCrackLibrary::DatasetIdOf("gap:1:2") == "gap");   // up to the first ':'
    }},
    {"unreadable_root", [] {
      BBRConfigManager::SetDataDir("/nonexistent/bbrsim-root");
      ExpectG4Exception("BBR011", [] { Pick("any", 500); }, "BBRCrackLibrary");
    }},
    {"concurrent_first_load", [] {
      TempDir d; Grid(d.path(), "cc", {"50", "500", "5000"});
      BBRConfigManager::SetDataDir(d.path().string());
      Pick("cc", 500);   // resolve the data root on the main thread first (std::thread is not a Geant4 worker)
      std::vector<std::thread> ts; std::vector<const void*> seen(8, nullptr);
      for (int i = 0; i < 8; ++i) ts.emplace_back([i, &seen] {
        double c; const auto& h = BBRCrackLibrary::Instance().Lookup("cc", (i % 3 == 0) ? 50 : (i % 3 == 1) ? 500 : 5000, c);
        seen[i] = &h;
      });
      for (auto& t : ts) t.join();
      std::set<const void*> distinct(seen.begin(), seen.end());
      CHECK(distinct.size() == 3);
    }},
    {"sidecar_required", [] {
      // One frequency of two lacks its sidecar: discovery stops even for a photon served by the other.
      TempDir d; Grid(d.path(), "two", {"50", "500"});
      std::filesystem::remove(d.path() / "waveguides" / "two_500GHz.dataset.json");
      BBRConfigManager::SetDataDir(d.path().string());
      ExpectG4Exception("BBR024", [] { Pick("two", 50); }, "BBRDatasetSidecar",
                        "two_500GHz.dataset.json: cannot be opened");
    }},
    {"sidecar_blocks_agree", [] {
      TempDir d; Grid(d.path(), "mix", {"50", "500"});
      auto p = hfssfix::SidecarFrom("mix_500GHz", hfssfix::Mini500());
      p.yHalf = 4.0e-3;   // another exit section at one frequency
      bbrtest::WriteFile(d.path() / "waveguides" / "mix_500GHz.dataset.json", hfssfix::SidecarJson(p));
      BBRConfigManager::SetDataDir(d.path().string());
      // The exit section, the modes and the extent_mm l it determines are named as the differing blocks.
      ExpectG4Exception("BBR024", [] { Pick("mix", 50); }, "BBRCrackLibrary::Discover",
                        "disagree on their frequency-independent physics, in block(s) geometry, modes, cross_section: ");
    }},
    {"validate_placed_cracks", [] {
      // Two placements of one id ("gap:1", "gap:2") with different solids: each is checked.
      TempDir d; Grid(d.path(), "gap", {"500"});
      BBRConfigManager::SetDataDir(d.path().string());
      auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* wlv = new G4LogicalVolume(new G4Box("W", 50 * mm, 50 * mm, 50 * mm), gal, "W");
      new G4PVPlacement(nullptr, {}, wlv, "W", nullptr, false, 0);
      auto* ok = new G4LogicalVolume(new G4Box("gapA", 2 * mm, 4.5 * mm, 0.025 * mm), BBRMaterials::GetVacuumWG(), "gapA");
      auto* thin = new G4LogicalVolume(new G4Box("gapB", 2 * mm, 5 * mm, 0.020 * mm), BBRMaterials::GetVacuumWG(), "gapB");
      new G4PVPlacement(nullptr, {0, 0, -10 * mm}, ok, "gap:1", wlv, false, 0);
      new G4PVPlacement(nullptr, {0, 0, 10 * mm}, thin, "gap:2", wlv, false, 0);
      ExpectG4Exception("BBR025", [] { BBRCrackLibrary::Instance().ValidatePlacedCracks(); }, "BBRDatasetSidecar",
                        "does not equal crack volume gap:2");
    }},
    {"validate_placed_cracks_ok", [] {
      TempDir d; Grid(d.path(), "gap", {"50", "500"});
      BBRConfigManager::SetDataDir(d.path().string());
      auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* wlv = new G4LogicalVolume(new G4Box("W", 50 * mm, 50 * mm, 50 * mm), gal, "W");
      new G4PVPlacement(nullptr, {}, wlv, "W", nullptr, false, 0);
      auto* ok = new G4LogicalVolume(new G4Box("gapA", 2 * mm, 4.5 * mm, 0.025 * mm), BBRMaterials::GetVacuumWG(), "gapA");
      new G4PVPlacement(nullptr, {}, ok, "gap", wlv, false, 0);
      BBRCrackLibrary::Instance().ValidatePlacedCracks();
      // Once per process: the second call is a no-op, so a crack placed since is not checked.
      auto* thin = new G4LogicalVolume(new G4Box("gapB", 2 * mm, 5 * mm, 0.020 * mm), BBRMaterials::GetVacuumWG(), "gapB");
      new G4PVPlacement(nullptr, {0, 0, 10 * mm}, thin, "gap:2", wlv, false, 0);
      BBRCrackLibrary::Instance().ValidatePlacedCracks();
      CHECK_NEAR(Pick("gap", 500), 500, 0);
    }},
    {"placed_cracks_json", [] {
      // One record per validated placement, in visiting order, with every grid frequency's sidecar.
      TempDir d; Grid(d.path(), "gap", {"50", "500"});
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK(BBRCrackLibrary::Instance().PlacedCracksJson() == "[]");
      auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* wlv = new G4LogicalVolume(new G4Box("W", 50 * mm, 50 * mm, 50 * mm), gal, "W");
      new G4PVPlacement(nullptr, {}, wlv, "W", nullptr, false, 0);
      auto* a = new G4LogicalVolume(new G4Box("gapA", 2 * mm, 4.5 * mm, 0.025 * mm), BBRMaterials::GetVacuumWG(), "gapA");
      auto* b = new G4LogicalVolume(new G4Box("gapB", 2 * mm, 4.5 * mm, 0.025 * mm), BBRMaterials::GetVacuumWG(), "gapB");
      new G4PVPlacement(nullptr, {0, 0, -10 * mm}, a, "gap:1", wlv, false, 0);
      new G4PVPlacement(nullptr, {0, 0, 10 * mm}, b, "gap:2", wlv, false, 0);
      CHECK(BBRCrackLibrary::Instance().PlacedCracksJson() == "[]");   // not validated yet
      BBRCrackLibrary::Instance().ValidatePlacedCracks();
      const auto arr = nlohmann::json::parse(BBRCrackLibrary::Instance().PlacedCracksJson());
      CHECK(arr.is_array() && arr.size() == 2);
      if (!arr.is_array() || arr.size() != 2) return;
      const std::vector<std::pair<std::string, std::vector<double>>> want = {
          {"gap:1", {4., 9., 0.05}}, {"gap:2", {4., 9., 0.05}}};   // told apart by name and placement
      for (std::size_t i = 0; i < 2; ++i) {
        const auto& c = arr[i];
        CHECK(c.at("volume") == want[i].first);
        CHECK(c.at("dataset_id") == "gap");
        CHECK(c.at("geant4_extent_mm").is_array() && c.at("geant4_extent_mm").size() == 3);
        for (std::size_t k = 0; k < 3; ++k) CHECK_NEAR(c.at("geant4_extent_mm").at(k).get<double>(), want[i].second[k], 1e-9);
        const auto& hx = c.at("hfss_extent_mm");
        CHECK_NEAR(hx.at("p").get<double>(), 4., 1e-12);
        CHECK_NEAR(hx.at("l").get<double>(), 9., 1e-12);
        CHECK_NEAR(hx.at("g").get<double>(), 0.05, 1e-12);
        const auto& fs = c.at("frequencies");
        CHECK(fs.is_array() && fs.size() == 2);
        if (!fs.is_array() || fs.size() != 2) continue;
        const char* labels[] = {"50GHz", "500GHz"};
        const double values[] = {50., 500.};
        for (std::size_t k = 0; k < 2; ++k) {
          CHECK(fs[k].at("label") == labels[k]);
          CHECK_NEAR(fs[k].at("frequency_ghz").get<double>(), values[k], 0);
          CHECK(fs[k].at("sidecar") == "gap_" + std::string(labels[k]) + ".dataset.json");
          CHECK(fs[k].at("recorded").is_object() && fs[k].at("recorded").at("frames").at("pose_rule") == "canonical-z");
        }
      }
    }},
    {"cutoff_warning", [] {
      TempDir d; Grid(d.path(), "one", {"500"}); Grid(d.path(), "low", {"10", "500"});
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK(Handler().CountWarnings("BBR026") == 0);
      CHECK_NEAR(Pick("one", 10), 500, 0);    // photon below cutoff, served above it: warns (a one-point grid too)
      CHECK_NEAR(Pick("one", 12), 500, 0);    // once per dataset and direction
      CHECK_NEAR(Pick("one", 600), 500, 0);   // both above: silent
      CHECK(Handler().CountWarnings("BBR026") == 1);
      CHECK_NEAR(Pick("low", 12), 10, 0);     // both below: silent
      CHECK_NEAR(Pick("low", 20), 10, 0);     // photon above, served by the 10 GHz point below: warns
      CHECK(Handler().CountWarnings("BBR026") == 2);
      CHECK(Handler().CountWarnings("BBR008") == 0);   // in-grid selections never clamp
      CHECK(Handler().CountWarnings("BBR026", "direction=below mode=TE10 ") == 1);
      CHECK(Handler().CountWarnings("BBR026", "direction=above mode=TE10 ") == 1);
    }},
    {"cutoff_warning_regime", [] {
      TempDir d; GridModes(d.path(), "c2", 500, false);
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK_NEAR(Pick("c2", 1600), 500, 0);   // above TE01 (1498.96), table below it; TE10 lies below both
      CHECK(Handler().CountWarnings("BBR026") == 1);
      CHECK(Handler().CountWarnings("BBR026", "dataset=c2 direction=above mode=TE01 ") == 1);
      CHECK_NEAR(Pick("c2", 1700), 500, 0);   // the same (mode, direction): silent
      CHECK(Handler().CountWarnings("BBR026") == 1);
      CHECK_NEAR(Pick("c2", 10), 500, 0);     // below TE10 (14.99), table above it; TE01 lies above both
      CHECK(Handler().CountWarnings("BBR026") == 2);
      CHECK(Handler().CountWarnings("BBR026", "dataset=c2 direction=below mode=TE10 ") == 1);
      CHECK_NEAR(Pick("c2", 10), 500, 0);
      CHECK(Handler().CountWarnings("BBR026") == 2);
    }},
    {"rho_limit_warning", [] {
      TempDir d;
      GridModes(d.path(), "disc2000", 2000, true);
      GridModes(d.path(), "disc1500", 1500, true);
      GridModes(d.path(), "disc1x", 2000, true, true);
      BBRConfigManager::SetDataDir(d.path().string());
      CHECK(Handler().CountWarnings("BBR028") == 0);
      CHECK_NEAR(Pick("disc1500", 1500), 1500, 0);   // below the 1756.98 limit
      CHECK_NEAR(Pick("disc1x", 2000), 2000, 0);     // schema 1: never
      CHECK(Handler().CountWarnings("BBR028") == 0);
      CHECK_NEAR(Pick("disc2000", 2000), 2000, 0);
      CHECK(Handler().CountWarnings("BBR028") == 1);
      CHECK(Handler().CountWarnings("BBR028", "BBR028 dataset=disc2000 grid_GHz=2000 valid_below_GHz=1756.98") == 1);
      CHECK_NEAR(Pick("disc2000", 2000), 2000, 0);
      CHECK(Handler().CountWarnings("BBR028") == 1);   // once per dataset
      CHECK(Handler().CountWarnings("BBR026") == 0);   // no regime cutoff lies between photon and table
    }},
  });
}
