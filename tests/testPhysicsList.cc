// testPhysicsList — BBSimPhysics in a physics list built as the examples build
// it (FTFP_BERT + EM option 4 + G4OpticalPhysics + BBSimPhysics), sequential run
// manager, trivial world, Initialize only (no event loop). Checks that the
// wrapper takes the stock boundary process's place in the optical photon's
// process list (D1) and that double or misordered registration is fatal
// (BBR014, D6).
#include "BBRTestSupport.hh"

#include "BBSimOpBoundaryProcess.hh"
#include "BBSimPhysics.hh"

#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalPhysics.hh"
#include "G4PVPlacement.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserDetectorConstruction.hh"

#include <string>
#include <vector>

using namespace bbrtest;

namespace {

class TrivialWorld : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* lv = new G4LogicalVolume(new G4Box("World", 1 * m, 1 * m, 1 * m), gal, "World");
    return new G4PVPlacement(nullptr, {}, lv, "World", nullptr, false, 0);
  }
};

// The optical photon's process list and post-step DoIt vector as names, with
// either boundary process (stock or wrapper) shown as "BOUNDARY".
struct Snapshot {
  std::vector<std::string> list, postDoIt;
  int nStock = 0, nWrapper = 0;   // top-level G4OpBoundaryProcess / BBSimOpBoundaryProcess in the list
};
std::string Token(const G4VProcess* p) {
  if (dynamic_cast<const BBSimOpBoundaryProcess*>(p) || dynamic_cast<const G4OpBoundaryProcess*>(p))
    return "BOUNDARY";
  return p->GetProcessName();
}
Snapshot Take() {
  Snapshot s;
  G4ProcessManager* pm = G4OpticalPhoton::Definition()->GetProcessManager();
  const G4ProcessVector* pl = pm->GetProcessList();
  for (std::size_t i = 0; i < pl->size(); ++i) {
    s.list.push_back(Token((*pl)[i]));
    if (dynamic_cast<const BBSimOpBoundaryProcess*>((*pl)[i])) ++s.nWrapper;
    else if (dynamic_cast<const G4OpBoundaryProcess*>((*pl)[i])) ++s.nStock;
  }
  const G4ProcessVector* pv = pm->GetPostStepProcessVector(typeDoIt);
  for (std::size_t i = 0; i < pv->size(); ++i) s.postDoIt.push_back(Token((*pv)[i]));
  return s;
}
std::string Join(const std::vector<std::string>& v) {
  std::string out;
  for (const auto& x : v) out += (out.empty() ? "" : " ") + x;
  return out;
}

// Registered between G4OpticalPhysics and BBSimPhysics: records the stock layout.
Snapshot gStock;
class Recorder : public G4VPhysicsConstructor {
 public:
  Recorder() : G4VPhysicsConstructor("BBRTestRecorder") {}
  void ConstructParticle() override {}
  void ConstructProcess() override { gStock = Take(); }
};

enum class Order { Normal, Twice, BBSimFirst };

G4RunManager* Build(Order order) {
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  rm->SetUserInitialization(new TrivialWorld);
  auto* pl = new FTFP_BERT(0);
  pl->ReplacePhysics(new G4EmStandardPhysics_option4());
  if (order == Order::BBSimFirst) pl->RegisterPhysics(new BBSimPhysics());
  pl->RegisterPhysics(new G4OpticalPhysics());
  if (order != Order::BBSimFirst) {
    pl->RegisterPhysics(new Recorder());
    pl->RegisterPhysics(new BBSimPhysics());
  }
  if (order == Order::Twice) pl->RegisterPhysics(new BBSimPhysics());
  rm->SetUserInitialization(pl);
  G4OpticalParameters::Instance()->SetProcessActivation("Cerenkov", false);
  G4OpticalParameters::Instance()->SetProcessActivation("Scintillation", false);
  return rm;
}

}  // namespace

int main(int argc, char** argv) {
  return RunCase(argc, argv, {
    {"wrapped_once_at_stock_index", [] {
      G4RunManager* rm = Build(Order::Normal);
      rm->Initialize();
      const Snapshot now = Take();
      std::printf("stock   list: %s\n        post: %s\n", Join(gStock.list).c_str(), Join(gStock.postDoIt).c_str());
      std::printf("wrapped list: %s\n        post: %s\n", Join(now.list).c_str(), Join(now.postDoIt).c_str());
      CHECK(gStock.nStock == 1 && gStock.nWrapper == 0);   // the recorder saw the stock layout
      CHECK(now.nStock == 0);                              // no top-level stock process left
      CHECK(now.nWrapper == 1);                            // exactly one wrapper
      CHECK(now.list == gStock.list);                      // same place in the process list
      CHECK(now.postDoIt == gStock.postDoIt);              // and in the post-step DoIt vector
      // The wrapper as Geant4 names it: G4WrapperProcess::RegisterProcess appends
      // the wrapped process's name, and the subtype keeps G4VProcess's default.
      const G4ProcessVector* pl = G4OpticalPhoton::Definition()->GetProcessManager()->GetProcessList();
      for (std::size_t i = 0; i < pl->size(); ++i) {
        if (!dynamic_cast<const BBSimOpBoundaryProcess*>((*pl)[i])) continue;
        CHECK((*pl)[i]->GetProcessName() == "BBSimOpBoundaryOpBoundary");
        CHECK((*pl)[i]->GetProcessSubType() == -1);
      }
    }},
    {"registered_twice_is_fatal", [] {
      G4RunManager* rm = Build(Order::Twice);
      ExpectG4Exception("BBR014", [&] { rm->Initialize(); }, "BBSimPhysics::WrapOpBoundaryProcess");
    }},
    {"wrong_order_is_fatal", [] {
      G4RunManager* rm = Build(Order::BBSimFirst);
      ExpectG4Exception("BBR014", [&] { rm->Initialize(); }, "BBSimPhysics::WrapOpBoundaryProcess");
    }},
  });
}
