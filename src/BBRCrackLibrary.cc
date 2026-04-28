#include "BBRCrackLibrary.hh"

BBRCrackLibrary& BBRCrackLibrary::Instance()
{
  static BBRCrackLibrary sInstance;
  return sInstance;
}

void BBRCrackLibrary::SetDataDir(const G4String& dir)
{
  fDataDir = dir;
}

BBRHFSSData& BBRCrackLibrary::Lookup(const G4String& datasetId)
{
  auto it = fCache.find(datasetId);
  if (it == fCache.end())
    it = fCache.emplace(datasetId,
                        std::make_unique<BBRHFSSData>(fDataDir, datasetId)).first;
  return *it->second;
}
