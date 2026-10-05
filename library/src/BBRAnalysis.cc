#include "BBRAnalysis.hh"

#include "BBRConfigManager.hh"
#include "BBRCrackLibrary.hh"
#include "G4AnalysisManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4TessellatedSolid.hh"
#include "Randomize.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Run.hh"
#include "G4MTRunManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdint>
#include "G4Version.hh"
#include "BBRBuildInfo.hh"

namespace {
// Default output file. It is registered with G4AnalysisManager::SetFileName in
// the ctor and opened by name-less OpenFile() at every run, so a macro can
// redirect a run with /analysis/setFileName (e.g. one file per run in a
// multi-run session; the same name is otherwise overwritten by the next run).
std::string JsonString(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') out << '\\' << c;
    else if (c < 0x20) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c);
    else out << c;
  }
  out << '"';
  return out.str();
}

const char* kOutputFile = "output/bbr.root";

// Distinct sentinels so an absent physical volume/material ("none", an expected
// case at the world boundary) is never confused with a name that was missing
// from the geometry stores at construction (kUnknownCode — a bug; it has no
// legend entry, so it surfaces as NaN in the Python loader rather than silently
// decoding to "none").
constexpr G4int kNoneCode = -1;
constexpr G4int kUnknownCode = -2;

// One table supplies the stable legend codes, stock status names, and event
// classification. Entries 0..20 retain their existing codes; new stock
// statuses are appended. -1 marks labels without a stock boundary enum.
struct StatusInfo {
  G4int stock;
  const char* name;
  G4int eventType;
};
#define STOCK_STATUS(status, type) {status, #status, type}
#define LOCAL_STATUS(name, type) {-1, name, type}
const std::vector<StatusInfo> kStatuses = {
    STOCK_STATUS(FresnelRefraction, 0), STOCK_STATUS(FresnelReflection, 1),
    {TotalInternalReflection, "TIR", 1}, STOCK_STATUS(LambertianReflection, 1),
    STOCK_STATUS(LobeReflection, 1), STOCK_STATUS(SpikeReflection, 1),
    STOCK_STATUS(BackScattering, 1), STOCK_STATUS(Absorption, 2),
    STOCK_STATUS(Detection, 2), STOCK_STATUS(NotAtBoundary, 3),
    STOCK_STATUS(SameMaterial, 3), STOCK_STATUS(StepTooSmall, 3),
    STOCK_STATUS(NoRINDEX, 2), LOCAL_STATUS("Other", 3),
    LOCAL_STATUS("BBRDiffractionTransmit", 0), LOCAL_STATUS("BBRDiffractionReflect", 1),
    LOCAL_STATUS("BBRReflect", 1), LOCAL_STATUS("BBRAbsorb", 2),
    LOCAL_STATUS("unknown", 3), LOCAL_STATUS("WorldExit", 3),
    LOCAL_STATUS("BulkAbsorption", 2),
    STOCK_STATUS(Undefined, 3), STOCK_STATUS(Transmission, 0),
    STOCK_STATUS(PolishedLumirrorAirReflection, 1),
    STOCK_STATUS(PolishedLumirrorGlueReflection, 1),
    STOCK_STATUS(PolishedAirReflection, 1),
    STOCK_STATUS(PolishedTeflonAirReflection, 1),
    STOCK_STATUS(PolishedTiOAirReflection, 1),
    STOCK_STATUS(PolishedTyvekAirReflection, 1),
    STOCK_STATUS(PolishedVM2000AirReflection, 1),
    STOCK_STATUS(PolishedVM2000GlueReflection, 1),
    STOCK_STATUS(EtchedLumirrorAirReflection, 1),
    STOCK_STATUS(EtchedLumirrorGlueReflection, 1),
    STOCK_STATUS(EtchedAirReflection, 1),
    STOCK_STATUS(EtchedTeflonAirReflection, 1),
    STOCK_STATUS(EtchedTiOAirReflection, 1),
    STOCK_STATUS(EtchedTyvekAirReflection, 1),
    STOCK_STATUS(EtchedVM2000AirReflection, 1),
    STOCK_STATUS(EtchedVM2000GlueReflection, 1),
    STOCK_STATUS(GroundLumirrorAirReflection, 1),
    STOCK_STATUS(GroundLumirrorGlueReflection, 1),
    STOCK_STATUS(GroundAirReflection, 1),
    STOCK_STATUS(GroundTeflonAirReflection, 1),
    STOCK_STATUS(GroundTiOAirReflection, 1),
    STOCK_STATUS(GroundTyvekAirReflection, 1),
    STOCK_STATUS(GroundVM2000AirReflection, 1),
    STOCK_STATUS(GroundVM2000GlueReflection, 1),
    STOCK_STATUS(Dichroic, 0),
    STOCK_STATUS(CoatedDielectricReflection, 1),
    STOCK_STATUS(CoatedDielectricRefraction, 0),
    STOCK_STATUS(CoatedDielectricFrustratedTransmission, 0)};
#undef STOCK_STATUS
#undef LOCAL_STATUS
}  // namespace

// FNV-1a over the name and bytes of every HFSS dataset sidecar
// (<data root>/waveguides/*.dataset.json). Each sidecar records the sha256 of
// its four CSVs (checked by validation/check_dataset_sidecars.py), so this
// identifies the data at a few kilobytes of reads instead of every CSV.
// "unavailable": no readable waveguides/ directory; "no-sidecars": one with
// no sidecar in it (the empty-input hash would look like a real fingerprint).
std::string BBRAnalysis::DataFingerprint(const std::string& dataRoot) {
  std::error_code ec;
  const auto root = std::filesystem::path(dataRoot) / "waveguides";
  if (!std::filesystem::is_directory(root, ec)) return "unavailable";
  const std::string suffix = ".dataset.json";
  std::vector<std::filesystem::path> paths;
  for (std::filesystem::directory_iterator it(root, ec), end; !ec && it != end; it.increment(ec)) {
    const std::string name = it->path().filename().string();
    if (it->is_regular_file() && name.size() > suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
      paths.push_back(it->path());
  }
  if (ec) return "unavailable";
  if (paths.empty()) return "no-sidecars";
  std::sort(paths.begin(), paths.end());
  std::uint64_t hash = 14695981039346656037ULL;
  auto add = [&](unsigned char c) { hash = (hash ^ c) * 1099511628211ULL; };
  for (const auto& p : paths) {
    for (unsigned char c : p.filename().string()) add(c);
    add(0);
    std::ifstream input(p, std::ios::binary);
    char buffer[65536];
    while (input.read(buffer, sizeof(buffer)) || input.gcount())
      for (std::streamsize i = 0; i < input.gcount(); ++i) add(static_cast<unsigned char>(buffer[i]));
    if (!input.eof()) return "unavailable";
    add(0);
  }
  std::ostringstream value;
  value << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << hash;
  return value.str();
}

BBRAnalysis::BBRAnalysis(const G4String& applicationVersion, const G4String& applicationFingerprint)
  : fApplicationVersion(applicationVersion), fApplicationFingerprint(applicationFingerprint) {
  BuildCategoryCodes();
  DefineNtuples();
}

BBRAnalysis::~BBRAnalysis() = default;

void BBRAnalysis::BuildCategoryCodes() {
  for (std::size_t i = 0; i < kStatuses.size(); ++i)
    fStatusCodes[kStatuses[i].name] = static_cast<G4int>(i);

  // Deterministic volume/material codes, sorted by name. Enumerate the
  // PHYSICAL volume store: the stepping/termination code records physical-volume
  // names (GetPhysicalVolume()->GetName()), so the legend must use the same
  // names or every vol_* / term_vol code decodes to kUnknownCode.
  std::vector<G4String> vols;
  for (const auto* pv : *G4PhysicalVolumeStore::GetInstance())
    vols.push_back(pv->GetName());
  std::sort(vols.begin(), vols.end());
  vols.erase(std::unique(vols.begin(), vols.end()), vols.end());
  for (std::size_t i = 0; i < vols.size(); ++i)
    fVolumeCodes[vols[i]] = static_cast<G4int>(i);
  fVolumeCodes["none"] = kNoneCode;

  std::vector<G4String> mats;
  for (const auto* m : *G4Material::GetMaterialTable())
    mats.push_back(m->GetName());
  std::sort(mats.begin(), mats.end());
  mats.erase(std::unique(mats.begin(), mats.end()), mats.end());
  for (std::size_t i = 0; i < mats.size(); ++i)
    fMaterialCodes[mats[i]] = static_cast<G4int>(i);
  fMaterialCodes["none"] = kNoneCode;
}

void BBRAnalysis::DefineNtuples() {
  auto* am = G4AnalysisManager::Instance();
  am->SetDefaultFileType("root");
  am->SetNtupleMerging(true);   // MT: merge worker ntuples into the master file
  am->SetVerboseLevel(1);

  // --- crossings ---
  am->SetFileName(kOutputFile);
  fCrossingsId = am->CreateNtuple("crossings", "Optical-photon boundary crossings");
  fCross.run_id = am->CreateNtupleIColumn("run_id");
  fCross.event_id = am->CreateNtupleIColumn("event_id");
  fCross.track_id = am->CreateNtupleIColumn("track_id");
  fCross.n_boundary = am->CreateNtupleIColumn("n_boundary");
  fCross.n_reflections = am->CreateNtupleIColumn("n_reflections");
  fCross.x = am->CreateNtupleDColumn("x_mm");
  fCross.y = am->CreateNtupleDColumn("y_mm");
  fCross.z = am->CreateNtupleDColumn("z_mm");
  fCross.energy = am->CreateNtupleDColumn("energy_eV");
  fCross.px_pre = am->CreateNtupleDColumn("px_pre");
  fCross.py_pre = am->CreateNtupleDColumn("py_pre");
  fCross.pz_pre = am->CreateNtupleDColumn("pz_pre");
  fCross.px_post = am->CreateNtupleDColumn("px_post");
  fCross.py_post = am->CreateNtupleDColumn("py_post");
  fCross.pz_post = am->CreateNtupleDColumn("pz_post");
  fCross.theta_in = am->CreateNtupleDColumn("theta_in_deg");
  fCross.phi_in = am->CreateNtupleDColumn("phi_in_deg");
  fCross.vol_pre = am->CreateNtupleIColumn("vol_pre_code");
  fCross.mat_pre = am->CreateNtupleIColumn("mat_pre_code");
  fCross.vol_post = am->CreateNtupleIColumn("vol_post_code");
  fCross.mat_post = am->CreateNtupleIColumn("mat_post_code");
  fCross.status = am->CreateNtupleIColumn("status_code");
  fCross.event_type = am->CreateNtupleIColumn("event_type_code");
  fCross.n_reflect = am->CreateNtupleIColumn("n_reflect");
  // HFSS grid frequency selected for a crack entry; -1 on every other crossing.
  fCross.hfss_freq = am->CreateNtupleDColumn("hfss_freq_GHz");
  am->FinishNtuple(fCrossingsId);

  // --- abspoints ---
  fAbsPointsId = am->CreateNtuple("abspoints", "Photon termination points");
  fAbs.run_id = am->CreateNtupleIColumn("run_id");
  fAbs.event_id = am->CreateNtupleIColumn("event_id");
  fAbs.track_id = am->CreateNtupleIColumn("track_id");
  fAbs.n_boundary = am->CreateNtupleIColumn("n_boundary");
  fAbs.n_reflections = am->CreateNtupleIColumn("n_reflections");
  fAbs.x = am->CreateNtupleDColumn("x_mm");
  fAbs.y = am->CreateNtupleDColumn("y_mm");
  fAbs.z = am->CreateNtupleDColumn("z_mm");
  fAbs.energy = am->CreateNtupleDColumn("energy_eV");
  fAbs.px = am->CreateNtupleDColumn("px");
  fAbs.py = am->CreateNtupleDColumn("py");
  fAbs.pz = am->CreateNtupleDColumn("pz");
  fAbs.n_reflect = am->CreateNtupleIColumn("n_reflect");
  fAbs.term_vol = am->CreateNtupleIColumn("term_vol_code");
  fAbs.term_status = am->CreateNtupleIColumn("term_status_code");
  am->FinishNtuple(fAbsPointsId);
}

void BBRAnalysis::BeginRun(const G4Run* run, G4bool master) {
  auto* am = G4AnalysisManager::Instance();
  std::filesystem::path root(std::string(am->GetFileName()));
  if (root.extension().empty()) root += ".root";
  if (root.extension() != ".root")
    G4Exception("BBRAnalysis::BeginRun", "BBR022", FatalException, "Output must be a ROOT file (.root)");
  fRootPath = root.string();
  std::error_code ec;
  if (!root.parent_path().empty()) std::filesystem::create_directories(root.parent_path(), ec);
  if (ec) G4Exception("BBRAnalysis::BeginRun", "BBR022", FatalException, ec.message().c_str());
  if (master) {
    auto metadata = root; metadata.replace_extension(".metadata.json");
    std::filesystem::remove(metadata, ec);
    if (ec) G4Exception("BBRAnalysis::BeginRun", "BBR022", FatalException, ec.message().c_str());
    BuildCategoryCodes();
    std::ostringstream config;
    config << std::setprecision(17);
    BBRConfigManager::Print(config);
    fConfiguration = config.str();
    // Geometry and the data directory are PreInit-only, so both are recorded once
    // per process: hashing the HFSS tables at every run would re-read them all.
    if (fGeometry.empty()) {
      std::ostringstream geometry;
      geometry << std::setprecision(17);
      for (const auto* pv : *G4PhysicalVolumeStore::GetInstance()) {
        const auto* material = pv->GetLogicalVolume()->GetMaterial();
        geometry << pv->GetName() << " copy=" << pv->GetCopyNo()
                 << " mother=" << (pv->GetMotherLogical() ? pv->GetMotherLogical()->GetName() : "none")
                 << " translation=" << pv->GetTranslation()
                 << " rotation=" << pv->GetObjectRotationValue()
                 << " material=" << (material ? material->GetName() : "none") << '\n';
        const auto* solid = pv->GetLogicalVolume()->GetSolid();
        // A tessellated (CAD) solid streams every facet; record its size instead.
        if (const auto* mesh = dynamic_cast<const G4TessellatedSolid*>(solid))
          geometry << "solid " << solid->GetName() << " type=" << solid->GetEntityType()
                   << " facets=" << mesh->GetNumberOfFacets() << '\n';
        else
          solid->StreamInfo(geometry);
      }
      fGeometry = geometry.str();
    }
    // The placed cracks and their sidecars: validated before the first run begins
    // (at /run/initialize, or in a sequential run manager's first BeamOn).
    if (fHfssDatasets.empty()) fHfssDatasets = BBRCrackLibrary::Instance().PlacedCracksJson();
    if (fDataFingerprint.empty())
      fDataFingerprint = DataFingerprint(BBRConfigManager::GetDataDir());
    std::ostringstream random;
    G4Random::getTheEngine()->put(random);
    fRandomState = random.str();
    BBRConfigManager::Print(G4cout);
  }
  if (!am->OpenFile())
    G4Exception("BBRAnalysis::BeginRun", "BBR022", FatalException, "Cannot open ROOT output");
  G4cout << "=== BBR Run " << run->GetRunID() << " begin (ROOT) ===" << G4endl;
}

void BBRAnalysis::WriteMetadata(const G4Run* run) const {
  auto path = std::filesystem::path(std::string(fRootPath));
  path.replace_extension(".metadata.json");
  auto temporary = path; temporary += ".tmp";
  std::ofstream js(temporary);
  auto dumpMap = [&](const std::map<G4String, G4int>& values) {
    js << "{";
    bool first = true;
    for (const auto& kv : values) {
      js << (first ? "" : ", ") << JsonString(std::to_string(kv.second)) << ": " << JsonString(kv.first);
      first = false;
    }
    js << "}";
  };
  const auto* mt = G4MTRunManager::GetMasterRunManager();
  js << "{\n  \"schema_version\": 2,\n  \"run_id\": " << run->GetRunID()
     << ",\n  \"events\": " << run->GetNumberOfEvent()
     << ",\n  \"worker_count\": " << (mt ? mt->GetNumberOfThreads() : 1)
     << ",\n  \"configuration\": " << JsonString(fConfiguration)
     << ",\n  \"geometry\": " << JsonString(fGeometry)
     << ",\n  \"random_engine_state\": " << JsonString(fRandomState)
     << ",\n  \"build\": {\"version\": " << JsonString(BBRSIM_BUILD_VERSION)
     << ", \"source_fingerprint\": " << JsonString(BBRSIM_SOURCE_FINGERPRINT)
     << ", \"application_version\": " << JsonString(fApplicationVersion)
     << ", \"application_source_fingerprint\": " << JsonString(fApplicationFingerprint)
     << ", \"geant4\": " << JsonString(G4Version) << "},\n  \"data\": {\"directory\": "
     << JsonString(std::filesystem::absolute(std::string(BBRConfigManager::GetDataDir())).string())
     << ", \"fingerprint\": " << JsonString(fDataFingerprint)
     << ", \"hfss_datasets\": " << fHfssDatasets
     << "},\n  \"legend\": {\n  \"status\": ";
  dumpMap(fStatusCodes);
  js << ",\n  \"event_type\": {\"0\": \"transmission\", \"1\": \"reflection\", \"2\": \"absorption\", \"3\": \"other\"},\n  \"volume\": ";
  dumpMap(fVolumeCodes);
  js << ",\n  \"material\": ";
  dumpMap(fMaterialCodes);
  js << "\n}}\n";
  js.close();
  if (!js) G4Exception("BBRAnalysis::WriteMetadata", "BBR022", FatalException, "Cannot write result metadata");
  std::error_code ec;
  std::filesystem::rename(temporary, path, ec);
  if (ec) G4Exception("BBRAnalysis::WriteMetadata", "BBR022", FatalException, ec.message().c_str());
}

void BBRAnalysis::EndRun(const G4Run* run, G4bool master) {
  auto* am = G4AnalysisManager::Instance();
  const auto written = am->Write();
  const auto closed = am->CloseFile();
  if (!written || !closed)
    G4Exception("BBRAnalysis::EndRun", "BBR022", FatalException, "Cannot write/close ROOT output");
  if (master) {
    WriteMetadata(run);
    G4cout << "=== BBR Run " << run->GetRunID() << " end: " << fRootPath << " written ===" << G4endl;
  }
}

G4int BBRAnalysis::EncodeStatus(const G4String& name) const {
  auto it = fStatusCodes.find(name);
  return it == fStatusCodes.end() ? fStatusCodes.at("unknown") : it->second;
}

G4int BBRAnalysis::EncodeVolume(const G4String& name) const {
  auto it = fVolumeCodes.find(name);
  return it == fVolumeCodes.end() ? kUnknownCode : it->second;
}

G4int BBRAnalysis::EncodeMaterial(const G4String& name) const {
  auto it = fMaterialCodes.find(name);
  return it == fMaterialCodes.end() ? kUnknownCode : it->second;
}

G4String BBRAnalysis::BoundaryStatusName(G4OpBoundaryProcessStatus status) {
  for (const auto& info : kStatuses)
    if (info.stock == status) return info.name;
  return "Other";
}

G4int BBRAnalysis::EventTypeForStatus(const G4String& name) {
  for (const auto& info : kStatuses)
    if (name == info.name) return info.eventType;
  return 3;
}
