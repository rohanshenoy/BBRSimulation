#ifndef BBRAnalysis_hh
#define BBRAnalysis_hh

#include "G4Run.hh"
#include "G4OpBoundaryProcess.hh"
#include "globals.hh"
#include <map>
#include <string>

// Owns the G4Analysis ROOT output (output/bbr.root) for a run.
// Ntuples: "crossings" (one row per optical-photon boundary crossing) and
// "abspoints" (one row per photon termination). The code->name dictionary is
// written by the master as <ROOT stem>.metadata.json (filling an
// ntuple on the master under SetNtupleMerging is unreliable).
// This helper is not a user action; each application retains its own actions.
// MT-safe: G4AnalysisManager keeps per-thread ntuples and merges them into the
// master file (SetNtupleMerging(true)). Categorical fields are stored as
// integer codes; codes are deterministic (built from the geometry at
// construction), so every worker thread computes the identical map with no lock.
class BBRAnalysis {
 public:
  explicit BBRAnalysis(const G4String& applicationVersion = "unspecified",
                       const G4String& applicationFingerprint = "unspecified");
  ~BBRAnalysis();

  void BeginRun(const G4Run* run, G4bool master);
  void EndRun(const G4Run* run, G4bool master);

  // Categorical encoders (string -> stable integer code). Volume/material
  // return -1 for "none" and -2 (kUnknownCode) for a name absent from the
  // legend; status falls back to the "unknown" code.
  G4int EncodeStatus(const G4String& name) const;
  G4int EncodeVolume(const G4String& name) const;
  G4int EncodeMaterial(const G4String& name) const;
  // Event-type code from a boundary-status name: 0=transmission, 1=reflection,
  // 2=absorption, 3=other.
  static G4String BoundaryStatusName(G4OpBoundaryProcessStatus status);
  static G4int EventTypeForStatus(const G4String& statusName);
  // Fingerprint of the HFSS data under <dataRoot>/waveguides: FNV-1a over the
  // top-level *.dataset.json sidecars only ("fnv1a64:" + 16 hex digits), else
  // "no-sidecars" or "unavailable" (no readable waveguides/ directory).
  static std::string DataFingerprint(const std::string& dataRoot);

  // Ntuple ids (assigned in the ctor; identical across threads).
  G4int fCrossingsId = -1;
  G4int fAbsPointsId = -1;

  // crossings column ids.
  struct CrossCols {
    G4int run_id, event_id, track_id, n_boundary, n_reflections, x, y, z, energy, px_pre, py_pre, pz_pre, px_post,
        py_post, pz_post, theta_in, phi_in, vol_pre, mat_pre, vol_post, mat_post,
        status, event_type, n_reflect, hfss_freq;
  } fCross;

  // abspoints column ids.
  struct AbsCols {
    G4int run_id, event_id, track_id, n_boundary, n_reflections, x, y, z, energy, px, py, pz, n_reflect, term_vol,
        term_status;
  } fAbs;

 private:
  void DefineNtuples();         // ctor: create the two ntuples
  void BuildCategoryCodes();    // ctor: status + volume + material -> codes
  void WriteMetadata(const G4Run* run) const; // master only, after successful CloseFile

  G4String fApplicationVersion;
  G4String fApplicationFingerprint;
  G4String fRootPath;
  G4String fConfiguration;
  G4String fGeometry;
  G4String fRandomState;
  G4String fDataFingerprint;
  G4String fHfssDatasets;       // BBRCrackLibrary::PlacedCracksJson(), raw JSON
  std::map<G4String, G4int> fStatusCodes;
  std::map<G4String, G4int> fVolumeCodes;
  std::map<G4String, G4int> fMaterialCodes;
};

#endif
