#ifndef BBRRunAction_hh
#define BBRRunAction_hh

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <map>

// Owns the G4Analysis ROOT output (output/bbr.root) for a run.
// Ntuples: "crossings" (one row per optical-photon boundary crossing) and
// "abspoints" (one row per photon termination). The code->name dictionary is
// written by the master as output/bbr_legend.json (a sidecar — filling an
// ntuple on the master under SetNtupleMerging is unreliable).
// MT-safe: G4AnalysisManager keeps per-thread ntuples and merges them into the
// master file (SetNtupleMerging(true)). Categorical fields are stored as
// integer codes; codes are deterministic (built from the geometry at
// construction), so every worker thread computes the identical map with no lock.
class BBRRunAction : public G4UserRunAction {
 public:
  BBRRunAction();
  ~BBRRunAction() override;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

  // Categorical encoders (string -> stable integer code, -1 if unknown).
  G4int EncodeStatus(const G4String& name) const;
  G4int EncodeVolume(const G4String& name) const;
  G4int EncodeMaterial(const G4String& name) const;
  // Event-type code from a boundary-status name: 0=transmission, 1=reflection,
  // 2=absorption, 3=other.
  static G4int EventTypeForStatus(const G4String& statusName);

  // Ntuple ids (assigned in the ctor; identical across threads).
  G4int fCrossingsId = -1;
  G4int fAbsPointsId = -1;

  // crossings column ids.
  struct CrossCols {
    G4int run_id, event_id, x, y, z, energy, px_pre, py_pre, pz_pre, px_post,
        py_post, pz_post, theta_in, phi_in, vol_pre, mat_pre, vol_post, mat_post,
        status, event_type, n_reflect;
  } fCross;

  // abspoints column ids.
  struct AbsCols {
    G4int run_id, event_id, x, y, z, energy, px, py, pz, n_reflect, term_vol,
        term_status;
  } fAbs;

 private:
  void DefineNtuples();         // ctor: create the two ntuples
  void BuildCategoryCodes();    // ctor: status + volume + material -> codes
  void WriteLegendJson() const; // master only, in BeginOfRun

  std::map<G4String, G4int> fStatusCodes;
  std::map<G4String, G4int> fVolumeCodes;
  std::map<G4String, G4int> fMaterialCodes;
};

#endif
