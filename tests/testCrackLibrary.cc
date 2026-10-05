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
#include <cmath>
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
      // The exit section, and the modes it determines, are named as the differing blocks.
      ExpectG4Exception("BBR024", [] { Pick("mix", 50); }, "BBRCrackLibrary::Discover",
                        "disagree on their frequency-independent physics, in block(s) modes, cross_section: ");
    }},
    {"validate_placed_cracks", [] {
      // Two placements of one id ("gap:1", "gap:2") with different solids: each is checked.
      TempDir d; Grid(d.path(), "gap", {"500"});
      BBRConfigManager::SetDataDir(d.path().string());
      auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* wlv = new G4LogicalVolume(new G4Box("W", 50 * mm, 50 * mm, 50 * mm), gal, "W");
      new G4PVPlacement(nullptr, {}, wlv, "W", nullptr, false, 0);
      auto* ok = new G4LogicalVolume(new G4Box("gapA", 2 * mm, 5 * mm, 0.026 * mm), BBRMaterials::GetVacuumWG(), "gapA");
      auto* thin = new G4LogicalVolume(new G4Box("gapB", 2 * mm, 5 * mm, 0.020 * mm), BBRMaterials::GetVacuumWG(), "gapB");
      new G4PVPlacement(nullptr, {0, 0, -10 * mm}, ok, "gap:1", wlv, false, 0);
      new G4PVPlacement(nullptr, {0, 0, 10 * mm}, thin, "gap:2", wlv, false, 0);
      ExpectG4Exception("BBR025", [] { BBRCrackLibrary::Instance().ValidatePlacedCracks(); }, "BBRDatasetSidecar",
                        "does not fit strictly inside crack volume gap:2");
    }},
    {"validate_placed_cracks_ok", [] {
      TempDir d; Grid(d.path(), "gap", {"50", "500"});
      BBRConfigManager::SetDataDir(d.path().string());
      auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* wlv = new G4LogicalVolume(new G4Box("W", 50 * mm, 50 * mm, 50 * mm), gal, "W");
      new G4PVPlacement(nullptr, {}, wlv, "W", nullptr, false, 0);
      auto* ok = new G4LogicalVolume(new G4Box("gapA", 2 * mm, 5 * mm, 0.026 * mm), BBRMaterials::GetVacuumWG(), "gapA");
      new G4PVPlacement(nullptr, {}, ok, "gap", wlv, false, 0);
      BBRCrackLibrary::Instance().ValidatePlacedCracks();
      // Once per process: the second call is a no-op, so a crack placed since is not checked.
      auto* thin = new G4LogicalVolume(new G4Box("gapB", 2 * mm, 5 * mm, 0.020 * mm), BBRMaterials::GetVacuumWG(), "gapB");
      new G4PVPlacement(nullptr, {0, 0, 10 * mm}, thin, "gap:2", wlv, false, 0);
      BBRCrackLibrary::Instance().ValidatePlacedCracks();
      CHECK_NEAR(Pick("gap", 500), 500, 0);
    }},
  });
}
