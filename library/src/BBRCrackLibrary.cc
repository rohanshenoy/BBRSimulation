#include "BBRCrackLibrary.hh"
#include "BBRConfigManager.hh"

#include "G4AutoLock.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ios.hh"

#include "nlohmann/json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <system_error>
#include <utility>
#include <vector>

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

G4String BBRCrackLibrary::DatasetIdOf(const G4String& volumeName)
{
  return volumeName.substr(0, volumeName.find(':'));
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

    set.entries.push_back({f, prefix + token + "GHz", nullptr, nullptr});
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

  // Every frequency's sidecar (fail closed), then the frequency-independent blocks must agree.
  for (auto& e : set.entries)
    e.sidecar = std::make_unique<BBRDatasetSidecar>(
        BBRDatasetSidecar::Load(fWaveguidesDir, datasetId, e.dirStem, e.freq_GHz));
  for (std::size_t i = 1; i < set.entries.size(); ++i) {
    const std::vector<std::string> differ = set.entries[i].sidecar->InvariantDiff(*set.entries[0].sidecar);
    if (!differ.empty()) {
      G4ExceptionDescription ed;
      ed << "The sidecars of dataset " << datasetId << " disagree on their frequency-independent physics, "
         << "in block(s)";
      for (std::size_t b = 0; b < differ.size(); ++b) ed << (b ? ", " : " ") << differ[b];
      ed << ": " << set.entries[0].sidecar->path << " and " << set.entries[i].sidecar->path << ".";
      G4Exception("BBRCrackLibrary::Discover", "BBR024", FatalException, ed);
    }
  }

  // Every frequency's four CSVs must be present before the first event, so an
  // incomplete tree stops the run at startup rather than mid-run. Existence
  // only: the CSVs are read on first selection (BBRHFSSData, same codes).
  for (const auto& e : set.entries) {
    for (const char* ephi : {"_Ephi=0", "_Ephi=1"}) {
      for (const auto& [file, code] : {std::pair{"far_field.csv", "BBR001"}, std::pair{"waveguide.csv", "BBR002"}}) {
        const fs::path p = dir / (e.dirStem + ephi) / file;
        std::error_code ec3;
        if (fs::is_regular_file(p, ec3)) continue;
        G4ExceptionDescription ed;
        ed << "Missing HFSS file " << AbsPath(p) << " (no such regular file): every frequency of dataset "
           << datasetId << " needs far_field.csv and waveguide.csv in both " << e.dirStem << "_Ephi=0 and "
           << e.dirStem << "_Ephi=1.";
        G4Exception("BBRCrackLibrary::Discover", code, FatalException, ed);
      }
    }
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

  // The lowest-mode cutoff (modes.cutoff_ghz): below it a closed guide carries
  // no propagating mode, so a table from the other side is the wrong physics.
  const G4double fc = E[k].sidecar->cutoffGHz;
  if (fc > 0. && ((nu_GHz < fc) != (E[k].freq_GHz < fc))) {
    const G4bool below = nu_GHz < fc;
    G4bool& warned = below ? set.warnedBelowCutoff : set.warnedAboveCutoff;
    if (!warned) {
      warned = true;
      G4ExceptionDescription ed;
      ed << "BBR026 dataset=" << datasetId << " direction=" << (below ? "below" : "above")
         << " nu_GHz=" << nu_GHz << " grid_GHz=" << E[k].freq_GHz << " cutoff_GHz=" << fc
         << " (" << E[k].sidecar->lowestMode << ") : the photon and the HFSS dataset serving it lie on "
         << "opposite sides of the lowest-mode cutoff; the table is used unchanged. Reported once per "
         << "dataset and direction.";
      G4Exception("BBRCrackLibrary::Lookup", "BBR026", JustWarning, ed);
    }
  }

  if (!E[k].data)
    E[k].data = std::make_unique<BBRHFSSData>(fWaveguidesDir, E[k].dirStem, E[k].freq_GHz, E[k].sidecar.get());

  chosen_GHz = E[k].freq_GHz;
  return *E[k].data;
}

void BBRCrackLibrary::ValidatePlacedCracks()
{
  G4AutoLock lock(&cacheMutex);
  if (fValidated) return;
  fValidated = true;
  for (const G4VPhysicalVolume* pv : *G4PhysicalVolumeStore::GetInstance()) {
    const G4LogicalVolume* lv = pv->GetLogicalVolume();
    const G4Material* mat = lv ? lv->GetMaterial() : nullptr;
    if (!mat || mat->GetName() != "vacuum_wg") continue;
    const G4String name = pv->GetName();
    const G4String id = DatasetIdOf(name);
    FrequencySet& set = Discover(id);
    const G4VSolid& solid = *lv->GetSolid();
    for (const auto& e : set.entries) e.sidecar->CheckFitsSolid(solid, name);
    G4ThreeVector lo, hi;
    solid.BoundingLimits(lo, hi);
    const G4ThreeVector extent = (hi - lo) / mm;
    fPlaced.push_back({name, id, {extent.x(), extent.y(), extent.z()}});
    const auto& sc = *set.entries.front().sidecar;
    G4cout << "[BBR] crack " << name << ": HFSS (p, l, g) = (" << sc.extentP_mm << ", " << sc.extentL_mm << ", "
           << sc.extentG_mm << ") mm, Geant4 (x, y, z) = (" << extent.x() << ", " << extent.y() << ", "
           << extent.z() << ") mm; " << set.entries.size() << " sidecar(s) fit" << G4endl;
  }
}

std::string BBRCrackLibrary::PlacedCracksJson() const
{
  using json = nlohmann::json;
  G4AutoLock lock(&cacheMutex);
  json out = json::array();
  for (const auto& p : fPlaced) {
    const FrequencySet& set = fSets.at(p.datasetId);   // discovered by ValidatePlacedCracks
    const auto& first = *set.entries.front().sidecar;
    json frequencies = json::array();
    for (const auto& e : set.entries)
      frequencies.push_back({{"label", e.sidecar->frequencyLabel},
                             {"frequency_ghz", e.freq_GHz},
                             {"sidecar", e.dirStem + ".dataset.json"},
                             {"recorded", e.sidecar->recorded.empty() ? json() : json::parse(e.sidecar->recorded)}});
    out.push_back({{"volume", std::string(p.volume)},
                   {"dataset_id", std::string(p.datasetId)},
                   {"geant4_extent_mm", {p.extent_mm[0], p.extent_mm[1], p.extent_mm[2]}},
                   {"hfss_extent_mm", {{"p", first.extentP_mm}, {"l", first.extentL_mm}, {"g", first.extentG_mm}}},
                   {"frequencies", frequencies}});
  }
  return out.dump(-1, ' ', false, json::error_handler_t::replace);
}
