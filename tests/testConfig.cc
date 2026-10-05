// testConfig.cc — BBRConfigManager / BBRConfigMessenger in a bare main (no
// run manager): defaults, the command table (states, broadcast), quoted
// /bbr/dataDir (D3), unit-aware setT with mK (D4), setter guards (D8), the
// data-dir precedence (one process per case), worker-clone isolation, the
// worker-before-master guard BBR020 (D5) and the /bbr/config/print flush (D17).
#include "BBRTestSupport.hh"

#include "BBRConfigManager.hh"

#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UImanager.hh"
#include "G4coutDestination.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <future>
#include <limits>
#include <string>
#include <thread>

#ifndef BBRSIM_TEST_PREFIX
#error "BBRSIM_TEST_PREFIX must be defined by tests/CMakeLists.txt"
#endif

namespace {

G4int Apply(const std::string& cmd) { return G4UImanager::GetUIpointer()->ApplyCommand(cmd); }

// Collects everything sent to G4cout and G4cerr while installed.
struct Capture : public G4coutDestination {
  std::string out, err;
  G4int ReceiveG4cout(const G4String& s) override { out += s; return 0; }
  G4int ReceiveG4cerr(const G4String& s) override { err += s; return 0; }
};

struct Cmd {
  const char* path;
  const char* arg;    // a valid argument
  bool idle;          // also available in Idle
  bool broadcast;
};

// The /bbr/ commands this test knows about. New commands are not an error.
const Cmd kCommands[] = {
  {"/bbr/config/print", "", true, false},
  {"/bbr/thermal/setT", "5", true, true},
  {"/bbr/thermal/emitterCenter", "-50 0 0 mm", true, true},
  {"/bbr/thermal/emitterSize", "1 20 20 mm", true, true},
  {"/bbr/gun/mode", "false", true, true},
  {"/bbr/gun/posX", "-20", true, true},
  {"/bbr/gun/posY", "0", true, true},
  {"/bbr/gun/posZ", "0", true, true},
  {"/bbr/gun/dirX", "1", true, true},
  {"/bbr/gun/dirY", "0", true, true},
  {"/bbr/gun/dirZ", "0", true, true},
  {"/bbr/gun/energy_eV", "2.07e-3", true, true},
  {"/bbr/gun/pol", "0 0 0", true, true},
  {"/bbr/dataDir", "/tmp/x", false, false},
  {"/bbr/det/setCuMaterial", "OFHC_Cu", false, false},
  {"/bbr/det/setCuRRR", "100", false, false},
  {"/bbr/det/setCuStageT", "4 K", false, false},
};

}  // namespace

int main(int argc, char** argv) {
  using namespace bbrtest;
  return RunCase(argc, argv, {
    {"defaults", [] {
      BBRConfigManager::Instance();
      CHECK(BBRConfigManager::GetThermalT_K() == 4.);
      CHECK(BBRConfigManager::GetEmitterCenter_mm() == G4ThreeVector(-50., 0., 0.));
      CHECK(BBRConfigManager::GetEmitterSize_mm() == G4ThreeVector(1., 20., 20.));
      CHECK(!BBRConfigManager::GetGunMode());
      CHECK(BBRConfigManager::GetGunPosX_mm() == -20. && BBRConfigManager::GetGunPosY_mm() == 0. &&
            BBRConfigManager::GetGunPosZ_mm() == 0.);
      CHECK(BBRConfigManager::GetGunDirX() == 1. && BBRConfigManager::GetGunDirY() == 0. &&
            BBRConfigManager::GetGunDirZ() == 0.);
      CHECK(BBRConfigManager::GetGunEnergy_eV() == 2.07e-3);
      CHECK(BBRConfigManager::GetGunPol() == G4ThreeVector(0., 0., 0.));
      CHECK(BBRConfigManager::GetCuRRR() == 100);
      CHECK(BBRConfigManager::GetCuStageT_K() == 4.);
    }},
    {"command_table", [] {
      BBRConfigManager::Instance();
      G4UIcommandTree* tree = G4UImanager::GetUIpointer()->GetTree();
      for (const auto& c : kCommands) {
        G4UIcommand* cmd = tree->FindPath(c.path);
        CHECK(cmd != nullptr);
        if (!cmd) continue;
        const auto* states = cmd->GetStateList();
        const bool preInit = std::find(states->begin(), states->end(), G4State_PreInit) != states->end();
        Report(preInit, __FILE__, __LINE__, std::string(c.path) + " allowed in PreInit");
        Report(cmd->ToBeBroadcasted() == c.broadcast, __FILE__, __LINE__,
               std::string(c.path) + " broadcast == " + (c.broadcast ? "true" : "false"));
      }
      CHECK(G4StateManager::GetStateManager()->SetNewState(G4State_Idle));
      for (const auto& c : kCommands) {
        const std::string line = std::string(c.path) + (*c.arg ? " " : "") + c.arg;
        const G4int rc = Apply(line);
        Report(rc == (c.idle ? 0 : 200), __FILE__, __LINE__,
               line + " in Idle: rc " + std::to_string(rc) + ", expected " + (c.idle ? "0" : "200"));
      }
    }},
    {"datadir_quotes", [] {  // D3
      BBRConfigManager::Instance();
      CHECK(Apply("/bbr/dataDir \"/tmp/a b\"") == 0);
      CHECK(BBRConfigManager::GetDataDir() == "/tmp/a b");
      CHECK(Apply("/bbr/dataDir '/tmp/c d'") == 0);
      CHECK(BBRConfigManager::GetDataDir() == "/tmp/c d");
      CHECK(Apply("/bbr/dataDir /tmp/e f") == 0);
      CHECK(BBRConfigManager::GetDataDir() == "/tmp/e f");
    }},
    {"setT_units", [] {  // D4
      BBRConfigManager::Instance();
      CHECK(Apply("/bbr/thermal/setT 10.0") == 0);
      CHECK(BBRConfigManager::GetThermalT_K() == 10.);
      CHECK(Apply("/bbr/thermal/setT 20 K") == 0);
      CHECK(BBRConfigManager::GetThermalT_K() == 20.);
      CHECK(Apply("/bbr/thermal/setT 300 mK") == 0);
      CHECK_REL(BBRConfigManager::GetThermalT_K(), 0.3, 1e-12);
      CHECK(Apply("/bbr/thermal/setT 10 mm") == 501);
      CHECK_REL(BBRConfigManager::GetThermalT_K(), 0.3, 1e-12);
      CHECK(Apply("/bbr/det/setCuStageT 4000 mK") == 0);
      CHECK_REL(BBRConfigManager::GetCuStageT_K(), 4., 1e-12);
      CHECK(Apply("/bbr/det/setCuStageT 40 kelvin") == 0);
      CHECK(BBRConfigManager::GetCuStageT_K() == 40.);
    }},
    {"validation_paths", [] {  // D8 for the setters
      BBRConfigManager::Instance();
      Capture cap;
      G4iosSetDestination(&cap);
      // A rejected value keeps the old one and fails the command.
      const G4int rc = Apply("/bbr/det/setCuRRR 0");
      CHECK(rc != 0 && BBRConfigManager::GetCuRRR() == 100);
      CHECK(cap.err.find("[BBR]") != std::string::npos);
      cap.err.clear();
      CHECK(Apply("/bbr/det/setCuMaterial bogus") != 0);
      CHECK(BBRConfigManager::GetCuRRR() == 100 && BBRConfigManager::GetCuStageT_K() == 4.);
      CHECK(cap.err.find("[BBR]") != std::string::npos);
      const double nan = std::numeric_limits<double>::quiet_NaN();
      CHECK(!BBRConfigManager::SetThermalT_K(nan));
      CHECK(BBRConfigManager::GetThermalT_K() == 4.);
      CHECK(!BBRConfigManager::SetCuStageT_K(nan));
      CHECK(BBRConfigManager::GetCuStageT_K() == 4.);
      CHECK(!BBRConfigManager::SetGunEnergy_eV(nan));
      CHECK(!BBRConfigManager::SetGunEnergy_eV(-1.));
      CHECK(BBRConfigManager::GetGunEnergy_eV() == 2.07e-3);
      CHECK(!BBRConfigManager::SetEmitterSize_mm(G4ThreeVector(nan, 1., 1.)));
      CHECK(BBRConfigManager::GetEmitterSize_mm() == G4ThreeVector(1., 20., 20.));
      CHECK(BBRConfigManager::SetCuRRR(100) && BBRConfigManager::SetCuMaterial("OFHC_Cu"));
      G4iosSetDestination(nullptr);
    }},
    {"datadir_default", [] {  // BBRSIMDATA unset: the compiled-in default
      CHECK(BBRConfigManager::GetDataDir() == BBRSIM_TEST_PREFIX "/share/BBRsim/data");
      Apply("/bbr/dataDir /cmd/data");
      CHECK(BBRConfigManager::GetDataDir() == "/cmd/data");
    }},
    {"datadir_env", [] {  // BBRSIMDATA=/env/data (CMake), then the command wins
      CHECK(BBRConfigManager::GetDataDir() == "/env/data");
      Apply("/bbr/dataDir /cmd/data");
      CHECK(BBRConfigManager::GetDataDir() == "/cmd/data");
    }},
    {"datadir_env_late", [] {  // BBRSIMDATA is read once, by the master constructor
      CHECK(BBRConfigManager::GetDataDir() == "/env/data");
      setenv("BBRSIMDATA", "/late/data", 1);
      CHECK(BBRConfigManager::GetDataDir() == "/env/data");
    }},
    {"worker_clone_isolation", [] {
      BBRConfigManager* master = BBRConfigManager::Instance();
      BBRConfigManager::SetThermalT_K(7.);
      BBRConfigManager::SetCuRRR(250);
      BBRConfigManager::SetCuStageT_K(40.);
      BBRConfigManager::SetDataDir("/master/data");
      std::promise<void> cloned, masterChanged;
      std::thread w([&] {
        G4Threading::G4SetThreadId(1);
        BBRConfigManager* inst = BBRConfigManager::Instance();
        CHECK(G4Threading::IsWorkerThread() && inst != master);
        CHECK(BBRConfigManager::GetThermalT_K() == 7. && BBRConfigManager::GetCuRRR() == 250);
        CHECK(BBRConfigManager::GetCuStageT_K() == 40.);
        CHECK(BBRConfigManager::GetDataDir() == "/master/data");
        BBRConfigManager::SetThermalT_K(33.);
        cloned.set_value();
        masterChanged.get_future().wait();
        CHECK(BBRConfigManager::GetThermalT_K() == 33.);   // its own change
        CHECK(BBRConfigManager::GetCuRRR() == 250);         // not the master's later 300
      });
      cloned.get_future().wait();
      CHECK(BBRConfigManager::GetThermalT_K() == 7.);      // the worker's change did not leak
      BBRConfigManager::SetCuRRR(300);
      masterChanged.set_value();
      w.join();
      CHECK(BBRConfigManager::GetCuRRR() == 300);
    }},
    {"worker_before_master", [] {  // D5: BBR020 instead of a null dereference
      std::string code, origin;
      std::thread w([&] {
        G4Threading::G4SetThreadId(1);
        ThrowingExceptionHandler h;  // the handler is per thread
        try {
          BBRConfigManager::Instance();
        } catch (const G4ExceptionCaught& e) {
          code = e.code;
          origin = e.origin;
        }
      });
      w.join();
      Report(code == "BBR020" && origin.find("BBRConfigManager::Instance") != std::string::npos,
             __FILE__, __LINE__, "worker-first Instance(): got '" + code + "' from '" + origin + "'");
    }},
    {"print_flushes", [] {  // D17
      BBRConfigManager::Instance();
      Capture cap;
      G4iosSetDestination(&cap);
      Apply("/bbr/config/print");
      const std::string got = cap.out;  // no G4endl after the command
      G4iosSetDestination(nullptr);
      CHECK(!got.empty());
      for (const char* f : {"/bbr/thermal/setT", "/bbr/thermal/emitterCenter", "/bbr/thermal/emitterSize",
                            "/bbr/gun/mode", "/bbr/gun/pos[XYZ]", "/bbr/gun/dir[XYZ]", "/bbr/gun/energy_eV",
                            "/bbr/gun/pol", "/bbr/dataDir", "/bbr/det/setCuRRR", "/bbr/det/setCuStageT"})
        Report(got.find(f) != std::string::npos, __FILE__, __LINE__, std::string("print shows ") + f);
    }},
  });
}
