#ifndef BBRCrackLibrary_hh
#define BBRCrackLibrary_hh

#include "BBRDatasetSidecar.hh"
#include "BBRHFSSData.hh"
#include "G4String.hh"

#include <map>
#include <memory>
#include <string>
#include <vector>

// Singleton that owns HFSS dataset discovery, loading and caching.
//
// Routing key: (crack volume name stripped of any ":N" placement suffix,
// nearest frequency). Datasets live under
// BBRConfigManager::GetDataDir()/waveguides as <id>_<freq>GHz_Ephi={0,1};
// the frequency grid for an id is discovered from those directory names on the
// first Lookup or ValidatePlacedCracks, and each frequency's CSVs are loaded on
// first selection.
// Adding a crack still requires only placing a new vacuum_wg volume whose name
// matches a dataset directory — no code changes.
//
// Discovery also loads every frequency's sidecar <id>_<freq>GHz.dataset.json
// (F1-F11, F13; fatal BBR024 or BBR025); the frequencies must agree on the
// frequency-independent blocks (BBR024, naming the blocks that differ). It
// then checks that both Ephi directories of every frequency hold far_field.csv
// (else fatal BBR001) and waveguide.csv (else BBR002) as regular files; the
// files are only checked to exist there, and are read on first selection.
//
// The cache is shared mutable state with synchronized lazy initialization, not
// an immutable singleton: discovery, selection, loading and the one-time clamp
// warnings all happen under one G4Mutex. Any future mutation must stay inside
// that lock. Entries are never erased, so a reference handed out by Lookup
// stays valid for the life of the process.
class BBRCrackLibrary
{
 public:
  static BBRCrackLibrary& Instance();

  // The HFSS dataset ID of a crack volume: its physical volume name up to the
  // first ':' ("gap:1" -> "gap", "gap" -> "gap").
  static G4String DatasetIdOf(const G4String& volumeName);

  // Dataset for `datasetId` at the grid frequency nearest nu_GHz in log space;
  // chosen_GHz receives that frequency. A photon outside the grid uses the
  // nearest edge and triggers one BBR008 warning per (dataset, side); a grid
  // with a single frequency never warns, since there is nothing to choose.
  // Fatal BBR011 if the data directory is unreadable or holds no
  // <id>_<freq>GHz_Ephi=0 directory for this id.
  const BBRHFSSData& Lookup(const G4String& datasetId, G4double nu_GHz,
                            G4double& chosen_GHz);

  // Validates every vacuum_wg volume of the geometry before the first event:
  // discovers its dataset (sidecars and CSV presence included) and checks that
  // each sidecar's exit cross-section fits strictly inside the volume's solid
  // (F12, fatal BBR025). Prints each crack's HFSS and Geant4 extents (F11). Called once per
  // process by BBSimOpBoundaryProcess::BuildPhysicsTable; later calls return.
  void ValidatePlacedCracks();

  // JSON array, one object per crack ValidatePlacedCracks checked, in its order:
  // volume, dataset_id, geant4_extent_mm, hfss_extent_mm and each frequency's
  // label, value, sidecar file and recorded fields; "[]" before it ran.
  std::string PlacedCracksJson() const;

 private:
  BBRCrackLibrary() = default;

  struct FrequencyEntry {
    G4double freq_GHz;                              // parsed from dirStem
    std::string dirStem;                            // "<id>_<token>GHz", token kept verbatim
    std::unique_ptr<BBRDatasetSidecar> sidecar;     // parsed at discovery (F1-F11, F13)
    std::unique_ptr<BBRHFSSData> data;              // loaded on first selection
  };

  struct FrequencySet {
    std::vector<FrequencyEntry> entries;  // ascending in freq_GHz
    G4bool warnedLow = false;
    G4bool warnedHigh = false;
  };

  // Scans the waveguide directory for this id. The caller must hold the lock.
  FrequencySet& Discover(const G4String& datasetId);

  struct PlacedCrack {
    G4String volume, datasetId;
    G4double extent_mm[3];                // solid's BoundingLimits extent along local x, y, z
  };

  G4String fWaveguidesDir;                  // <data root>/waveguides, resolved once
  std::map<G4String, FrequencySet> fSets;
  G4bool fValidated = false;                // ValidatePlacedCracks ran
  std::vector<PlacedCrack> fPlaced;         // the cracks it checked, in its order
};

#endif
