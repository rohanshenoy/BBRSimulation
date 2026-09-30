#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"
#include "BBRConfigManager.hh"
#include "BBRCrackLibrary.hh"
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
      ExpectG4Exception("BBR001", [] { Pick("half", 500); }, "BBRHFSSData");        // Ephi=1 missing: at first selection
      ExpectG4Exception("BBR001", [] { Pick("half", 500); }, "BBRHFSSData");        // and again (the lock was released)
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
  });
}
