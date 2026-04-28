#ifndef BBRCrackLibrary_hh
#define BBRCrackLibrary_hh

#include "BBRHFSSData.hh"
#include "G4String.hh"

#include <map>
#include <memory>

// Singleton that owns HFSS dataset loading and caching.
// Routing key: volume name stripped of any ":N" placement suffix = dataset ID.
// Lazy-loads BBRHFSSData on first Lookup; safe for single-threaded runs.
class BBRCrackLibrary
{
 public:
  static BBRCrackLibrary& Instance();

  void         SetDataDir(const G4String& dir);    // default: "../HFSSSimData"
  BBRHFSSData& Lookup(const G4String& datasetId);  // lazy-loads on first call

 private:
  BBRCrackLibrary() = default;
  G4String fDataDir = "../HFSSSimData";
  std::map<G4String, std::unique_ptr<BBRHFSSData>> fCache;
};

#endif
