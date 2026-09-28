#include "BBRCrackLibrary.hh"
#include "BBRConfigManager.hh"

#include "G4AutoLock.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <system_error>

namespace {
G4Mutex cacheMutex = G4MUTEX_INITIALIZER;

// Absolute path for error messages, falling back to the input if the working
// directory cannot be resolved.
std::string AbsPath(const std::filesystem::path& p)
{
  std::error_code ec;
  const auto abs = std::filesystem::absolute(p, ec);
  return ec ? p.string() : abs.string();
}
}  // namespace

BBRCrackLibrary& BBRCrackLibrary::Instance()
{
  static BBRCrackLibrary sInstance;
  return sInstance;
}

BBRCrackLibrary::FrequencySet& BBRCrackLibrary::Discover(const G4String& datasetId)
{
  auto it = fSets.find(datasetId);
  if (it != fSets.end()) return it->second;

  namespace fs = std::filesystem;
  if (fWaveguidesDir.empty())
    fWaveguidesDir = BBRConfigManager::GetDataDir() + "/waveguides";
  const char* env = std::getenv("BBRSIMDATA");
  const fs::path dir(fWaveguidesDir.c_str());

  // The error_code overload: a missing or unreadable directory must surface as
  // a G4Exception, not as a C++ exception escaping a worker thread.
  std::error_code ec;
  fs::directory_iterator dit(dir, ec);
  if (ec) {
    G4ExceptionDescription ed;
    ed << "Cannot read the HFSS waveguide directory " << AbsPath(dir) << " ("
       << ec.message() << "). Set /bbr/dataDir <root> before /run/initialize, "
       << "where <root> contains waveguides/, or set BBRSIMDATA (currently "
       << (env ? env : "unset") << ").";
    G4Exception("BBRCrackLibrary::Discover", "BBR011", FatalException, ed);
  }

  const std::string prefix = std::string(datasetId) + "_";
  const std::string suffix = "GHz_Ephi=0";
  const std::string legacy = std::string(datasetId) + "_Ephi=0";

  FrequencySet set;
  for (const auto& entry : dit) {
    std::error_code ec2;
    if (!entry.is_directory(ec2) || ec2) continue;
    const std::string name = entry.path().filename().string();

    if (name == legacy) {
      G4ExceptionDescription ed;
      ed << "Dataset directory " << AbsPath(entry.path())
         << " carries no frequency. Rename it to " << prefix << "<freq>" << suffix
         << " (e.g. " << prefix << "500" << suffix << "); the lookup is keyed by "
         << "frequency.";
      G4Exception("BBRCrackLibrary::Discover", "BBR011", FatalException, ed);
    }
    if (name.size() <= prefix.size() + suffix.size()) continue;
    if (name.compare(0, prefix.size(), prefix) != 0) continue;
    if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;

    // The token must parse completely; "<id>_junkGHz_Ephi=0" is not a dataset.
    const std::string token =
        name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
    char* end = nullptr;
    const double f = std::strtod(token.c_str(), &end);
    if (token.empty() || end == nullptr || *end != '\0' || !(f > 0.)) continue;

    set.entries.push_back({f, prefix + token + "GHz", nullptr});
  }

  std::sort(set.entries.begin(), set.entries.end(),
            [](const FrequencyEntry& a, const FrequencyEntry& b) {
              return a.freq_GHz < b.freq_GHz;
            });

  // Distinct tokens can parse to the same value ("500" and "5e2"); that would
  // make the selected dataset depend on directory order, so it is an error.
  for (std::size_t i = 1; i < set.entries.size(); ++i) {
    if (set.entries[i].freq_GHz == set.entries[i - 1].freq_GHz) {
      G4ExceptionDescription ed;
      ed << "Duplicate frequency " << set.entries[i].freq_GHz << " GHz for dataset "
         << datasetId << ": directories " << set.entries[i - 1].dirStem << " and "
         << set.entries[i].dirStem << " name the same frequency.";
      G4Exception("BBRCrackLibrary::Discover", "BBR011", FatalException, ed);
    }
  }

  if (set.entries.empty()) {
    G4ExceptionDescription ed;
    ed << "No " << prefix << "<freq>" << suffix << " directory under " << AbsPath(dir)
       << " for crack volume " << datasetId
       << ". Set /bbr/dataDir <root> (root contains waveguides/) or BBRSIMDATA "
       << "(currently " << (env ? env : "unset") << ").";
    G4Exception("BBRCrackLibrary::Discover", "BBR011", FatalException, ed);
  }

  G4cout << "[BBR] HFSS dataset " << datasetId << ": " << set.entries.size()
         << " frequency grid point(s):";
  for (const auto& e : set.entries) G4cout << " " << e.freq_GHz;
  G4cout << " GHz (from " << AbsPath(dir) << ")" << G4endl;

  return fSets.emplace(datasetId, std::move(set)).first->second;
}

const BBRHFSSData& BBRCrackLibrary::Lookup(const G4String& datasetId, G4double nu_GHz,
                                           G4double& chosen_GHz)
{
  // The lock covers discovery, selection, the warning flags and the lazy load:
  // worker threads call this concurrently and all of it mutates shared state.
  G4AutoLock lock(&cacheMutex);
  FrequencySet& set = Discover(datasetId);
  auto& E = set.entries;

  // Selection rule: nearest grid point in log frequency, ties to the lower,
  // clamped at the edges. Keep in lock-step with bbrsim.hfss.select_frequency,
  // which the validators compare against.
  std::size_t k = 0;
  G4int clamped = 0;
  if (E.size() == 1) {
    k = 0;                                  // single-frequency grid: no choice, no warning
  } else if (nu_GHz < E.front().freq_GHz) {
    k = 0;
    clamped = -1;
  } else if (nu_GHz > E.back().freq_GHz) {
    k = E.size() - 1;
    clamped = +1;
  } else {
    const G4double lnu = std::log10(nu_GHz);
    G4double best = std::numeric_limits<G4double>::infinity();
    for (std::size_t i = 0; i < E.size(); ++i) {
      // Ascending order + strict '<' keeps the LOWER frequency on an exact tie.
      const G4double d = std::abs(std::log10(E[i].freq_GHz) - lnu);
      if (d < best) { best = d; k = i; }
    }
  }

  if (clamped != 0) {
    G4bool& warned = (clamped < 0) ? set.warnedLow : set.warnedHigh;
    if (!warned) {
      warned = true;
      G4ExceptionDescription ed;
      ed << "BBR008 dataset=" << datasetId << " side=" << (clamped < 0 ? "low" : "high")
         << " nu_GHz=" << nu_GHz << " edge_GHz=" << E[k].freq_GHz
         << " : photon frequency lies outside the HFSS grid; the edge dataset is "
         << "used unchanged. Reported once per dataset and side.";
      G4Exception("BBRCrackLibrary::Lookup", "BBR008", JustWarning, ed);
    }
  }

  if (!E[k].data)
    E[k].data = std::make_unique<BBRHFSSData>(fWaveguidesDir, E[k].dirStem, E[k].freq_GHz);

  chosen_GHz = E[k].freq_GHz;
  return *E[k].data;
}
