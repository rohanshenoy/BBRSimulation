// Consumer smoke test: compiles against the installed BBRsim headers, links
// libBBRsim through RPATH and prints the master data default.
#include "BBRConfigManager.hh"
#include "BBSimPhysics.hh"
#include "G4ios.hh"

int main() {
  G4cout << "BBRsim data dir: " << BBRConfigManager::GetDataDir() << G4endl;
  return 0;
}
