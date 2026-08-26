#include "BBRCrackLibrary.hh"

#include "G4AutoLock.hh"

namespace {
G4Mutex cacheMutex = G4MUTEX_INITIALIZER;
}

BBRCrackLibrary& BBRCrackLibrary::Instance()
{
  static BBRCrackLibrary sInstance;
  return sInstance;
}

void BBRCrackLibrary::SetDataDir(const G4String& dir)
{
  G4AutoLock lock(&cacheMutex);
  fDataDir = dir;
}

BBRHFSSData& BBRCrackLibrary::Lookup(const G4String& datasetId)
{
  // The cache is shared by worker threads. Hold the lock through construction
  // so no worker can read or mutate the map while another worker loads a CSV
  // dataset or while the data directory is being changed.
  G4AutoLock lock(&cacheMutex);
  auto it = fCache.find(datasetId);
  if (it == fCache.end())
    it = fCache.emplace(datasetId,
                        std::make_unique<BBRHFSSData>(fDataDir, datasetId)).first;
  return *it->second;
}
