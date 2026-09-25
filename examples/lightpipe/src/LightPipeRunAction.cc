#include "LightPipeRunAction.hh"

#include "BBRConfigManager.hh"
#include "G4AnalysisManager.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>

namespace {
const char* kOutputDir = "output";
// Default output file. It is registered with G4AnalysisManager::SetFileName in
// the ctor and opened by name-less OpenFile() at every run, so a macro can
// redirect a run with /analysis/setFileName (e.g. one file per run in a
// multi-run session; the same name is otherwise overwritten by the next run).
const char* kOutputFile = "output/bbr.root";
const char* kLegendFile = "output/bbr_legend.json";

// Distinct sentinels so an absent physical volume/material ("none", an expected
// case at the world boundary) is never confused with a name that was missing
// from the geometry stores at construction (kUnknownCode — a bug; it has no
// legend entry, so it surfaces as NaN in the Python loader rather than silently
// decoding to "none").
constexpr G4int kNoneCode = -1;
constexpr G4int kUnknownCode = -2;

// Stock + BBR-wrapper boundary statuses, in a fixed order so their codes are
// stable across builds (new names are appended, never inserted). Keep in sync
// with LightPipeSteppingAction::StatusStr / TerminationStatus and the
// BBSimOpBoundaryProcess status strings. "WorldExit" and "BulkAbsorption" are
// termination-only labels (abspoints.term_status).
const std::vector<G4String> kStatuses = {
    "FresnelRefraction", "FresnelReflection", "TIR", "LambertianReflection",
    "LobeReflection", "SpikeReflection", "BackScattering", "Absorption",
    "Detection", "NotAtBoundary", "SameMaterial", "StepTooSmall", "NoRINDEX",
    "Other", "BBRDiffractionTransmit", "BBRDiffractionReflect", "BBRReflect",
    "BBRAbsorb", "unknown", "WorldExit", "BulkAbsorption"};
}  // namespace

LightPipeRunAction::LightPipeRunAction() : G4UserRunAction() {
  BuildCategoryCodes();
  DefineNtuples();
}

LightPipeRunAction::~LightPipeRunAction() = default;

void LightPipeRunAction::BuildCategoryCodes() {
  for (std::size_t i = 0; i < kStatuses.size(); ++i)
    fStatusCodes[kStatuses[i]] = static_cast<G4int>(i);

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

void LightPipeRunAction::DefineNtuples() {
  auto* am = G4AnalysisManager::Instance();
  am->SetDefaultFileType("root");
  am->SetNtupleMerging(true);   // MT: merge worker ntuples into the master file
  am->SetVerboseLevel(1);

  // --- crossings ---
  am->SetFileName(kOutputFile);
  fCrossingsId = am->CreateNtuple("crossings", "Optical-photon boundary crossings");
  fCross.run_id = am->CreateNtupleIColumn("run_id");
  fCross.event_id = am->CreateNtupleIColumn("event_id");
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

void LightPipeRunAction::BeginOfRunAction(const G4Run* run) {
  std::error_code ec;
  std::filesystem::create_directories(kOutputDir, ec);

  auto* am = G4AnalysisManager::Instance();
  am->OpenFile();   // kOutputFile unless /analysis/setFileName changed it

  // Under MT, BuildForMaster() runs the master's ctor before the master builds
  // the detector, so the geometry stores were empty and the master's
  // volume/material maps held only "none". By BeginOfRunAction the master's
  // geometry is constructed, so rebuild the maps here before emitting the
  // legend. The enumeration is deterministic (sorted by name), so the master's
  // codes match the codes the workers fill from their own (identical) stores.
  if (IsMaster()) {
    BuildCategoryCodes();
    WriteLegendJson();
    BBRConfigManager::Print(G4cout);   // run provenance
  }
  G4cout << "=== BBR Run " << run->GetRunID() << " begin (ROOT) ===" << G4endl;
}

// Master-only. Writes {category: {code: name}}. Geant4 volume/material/status
// names are bare identifiers (no quotes/backslashes), so no JSON escaping is
// needed. A sidecar avoids the unreliable master-thread ntuple fill under
// SetNtupleMerging.
void LightPipeRunAction::WriteLegendJson() const {
  std::ofstream js(kLegendFile);
  auto dumpMap = [&](const std::map<G4String, G4int>& m) {
    js << "{";
    bool first = true;
    for (const auto& kv : m) {
      js << (first ? "" : ", ") << "\"" << kv.second << "\": \"" << kv.first
         << "\"";
      first = false;
    }
    js << "}";
  };
  js << "{\n  \"status\": ";
  dumpMap(fStatusCodes);
  js << ",\n  \"event_type\": {\"0\": \"transmission\", \"1\": \"reflection\", "
        "\"2\": \"absorption\", \"3\": \"other\"},\n  \"volume\": ";
  dumpMap(fVolumeCodes);
  js << ",\n  \"material\": ";
  dumpMap(fMaterialCodes);
  js << "\n}\n";
}

void LightPipeRunAction::EndOfRunAction(const G4Run* run) {
  auto* am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();
  if (IsMaster())
    G4cout << "=== BBR Run " << run->GetRunID() << " end: " << am->GetFileName()
           << " written ===" << G4endl;
}

G4int LightPipeRunAction::EncodeStatus(const G4String& name) const {
  auto it = fStatusCodes.find(name);
  return it == fStatusCodes.end() ? fStatusCodes.at("unknown") : it->second;
}

G4int LightPipeRunAction::EncodeVolume(const G4String& name) const {
  auto it = fVolumeCodes.find(name);
  return it == fVolumeCodes.end() ? kUnknownCode : it->second;
}

G4int LightPipeRunAction::EncodeMaterial(const G4String& name) const {
  auto it = fMaterialCodes.find(name);
  return it == fMaterialCodes.end() ? kUnknownCode : it->second;
}

G4int LightPipeRunAction::EventTypeForStatus(const G4String& s) {
  if (s == "FresnelRefraction" || s == "BBRDiffractionTransmit") return 0;
  if (s == "FresnelReflection" || s == "TIR" || s == "SpikeReflection" ||
      s == "LobeReflection" || s == "LambertianReflection" ||
      s == "BackScattering" || s == "BBRDiffractionReflect" || s == "BBRReflect")
    return 1;
  if (s == "Absorption" || s == "Detection" || s == "BBRAbsorb" ||
      s == "BulkAbsorption")
    return 2;
  return 3;  // NotAtBoundary, SameMaterial, StepTooSmall, NoRINDEX, Other,
             // unknown, WorldExit
}
