// testPhysicsList — BBSimPhysics in a physics list built as the examples build
// it (FTFP_BERT + EM option 4 + G4OpticalPhysics + BBSimPhysics). The first
// three cases use a sequential run manager and a trivial world and stop at
// Initialize (no event loop). They check that the wrapper takes the stock
// boundary process's place in the optical photon's process list (D1) and that
// double or misordered registration is fatal (BBR014, D6).
//
// The event_loop_fatal cases run one event (optical physics only) into a
// vacuum_wg crack whose dataset passes the startup check (sidecar and all four
// CSVs present) but whose Ephi=1 waveguide.csv has a non-numeric field in a
// data row: BBR013, raised by BBRHFSSData when it reads the CSVs at the first
// selection under the wrapper's PostStepDoIt, must leave BeamOn as the
// G4Exception. With the sequential run manager the throwing handler catches
// it; with the tasking run manager it is raised on a worker thread, whose own
// stock handler aborts the process (event_loop_fatal_mt, run through
// tests/ExpectAbort.cmake).
//
// The *_stops_before_events cases: BBSimOpBoundaryProcess::BuildPhysicsTable
// validates every placed crack at the first run initialization (/run/beamOn
// with the sequential run manager used there; /run/initialize with an MT or
// task run manager, whose Initialize calls BeamOn(0)), so a dataset without
// its sidecar (BBR024) or with an exit section wider than the crack (BBR025)
// stops the run before any event.
#include "BBRTestSupport.hh"
#include "HFSSFixture.hh"

#include "BBRConfigManager.hh"
#include "BBRMaterials.hh"
#include "BBSimOpBoundaryProcess.hh"
#include "BBSimPhysics.hh"

#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalPhysics.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <filesystem>
#include <string>
#include <utility>
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

// The event-loop world: a vacuum_wg crack (the HFSS dataset id is the volume
// name), entry face at x = 0, in a G4_Galactic world.
class CrackWorld : public G4VUserDetectorConstruction {
 public:
  explicit CrackWorld(std::string name) : fName(std::move(name)) {}
  G4VPhysicalVolume* Construct() override {
    auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* wlv = new G4LogicalVolume(new G4Box("World", 50 * mm, 50 * mm, 50 * mm), gal, "World");
    auto* world = new G4PVPlacement(nullptr, {}, wlv, "World", nullptr, false, 0);
    auto* clv = new G4LogicalVolume(new G4Box(fName, 2 * mm, 5 * mm, 0.026 * mm),
                                    BBRMaterials::GetVacuumWG(), fName);
    new G4PVPlacement(nullptr, {2 * mm, 0, 0}, clv, fName, wlv, false, 0);
    return world;
  }
 private:
  std::string fName;
};

// Data root for the event-loop cases: ./hfss in the case's working directory.
std::filesystem::path DataRoot() { return std::filesystem::current_path() / "hfss"; }

// A complete dataset (sidecar and all four CSVs) whose Ephi=1 waveguide.csv
// carries Ex_real = "abc" in its data row under a valid header: the startup
// check, which looks only for the files, passes, and the first selection,
// mid-event, raises BBR013 when the loader parses that row.
void WriteBadRowDataset(const std::string& id) {
  auto ds = hfssfix::Mini500();
  hfssfix::WriteDataset(DataRoot(), id + "_500GHz", ds);   // the sidecar is derived from the valid rows
  ds.wg1 = std::string(hfssfix::WG_HDR) + "500GHz,1,0,180,1,4,0,0,0,abc,0,1,0,0,0\n";
  bbrtest::WriteFile(DataRoot() / "waveguides" / (id + "_500GHz_Ephi=1") / "waveguide.csv", ds.wg1);
}

// One 500 GHz photon from x = -20 mm along +x, into the crack's entry face.
class CrackGun : public G4VUserPrimaryGeneratorAction {
 public:
  CrackGun() : fGun(1) {
    fGun.SetParticleDefinition(G4OpticalPhoton::Definition());
    fGun.SetParticleEnergy(2.07e-3 * eV);
    fGun.SetParticlePosition({-20 * mm, 0, 0});
    fGun.SetParticleMomentumDirection({1, 0, 0});
    fGun.SetParticlePolarization({0, 1, 0});
  }
  void GeneratePrimaries(G4Event* ev) override { fGun.GeneratePrimaryVertex(ev); }
 private:
  G4ParticleGun fGun;
};

class CrackActions : public G4VUserActionInitialization {
 public:
  void Build() const override { SetUserAction(new CrackGun); }
};

// Optical physics only: an event with the examples' FTFP_BERT list would first
// build every physics table, which needs most of the Geant4 datasets.
G4RunManager* BuildEventLoop(G4RunManagerType type, const std::string& crack) {
  BBRConfigManager::SetDataDir(DataRoot().string());   // on the master, before the workers copy it
  auto* rm = G4RunManagerFactory::CreateRunManager(type);
  bbrtest::CreateElectronInPreInit();
  rm->SetNumberOfThreads(2);            // no effect on the sequential run manager
  rm->SetUserInitialization(new CrackWorld(crack));
  G4OpticalParameters::Instance()->SetProcessActivation("Cerenkov", false);
  G4OpticalParameters::Instance()->SetProcessActivation("Scintillation", false);
  auto* pl = new G4VModularPhysicsList;
  pl->RegisterPhysics(new G4OpticalPhysics());
  pl->RegisterPhysics(new BBSimPhysics());
  rm->SetUserInitialization(pl);
  rm->SetUserInitialization(new CrackActions);
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
    {"event_loop_fatal", [] {
      WriteBadRowDataset("BadRowCrack");
      G4RunManager* rm = BuildEventLoop(G4RunManagerType::SerialOnly, "BadRowCrack");
      rm->Initialize();
      ExpectG4Exception("BBR013", [&] { rm->BeamOn(1); }, "BBRHFSSData::Load",
                        "Bad numeric field Ex_real = 'abc' in " + DataRoot().string() +
                          "/waveguides/BadRowCrack_500GHz_Ephi=1/waveguide.csv");
      // Raised mid-event, on the step that entered the crack. The run is left
      // half-processed, so this is the last thing the process does.
      CHECK(G4StateManager::GetStateManager()->GetCurrentState() == G4State_EventProc);
      const G4Step* s = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager()->GetStep();
      CHECK(s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "World");
      CHECK(s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "BadRowCrack");
      CHECK(s->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
      CHECK_NEAR(s->GetPostStepPoint()->GetPosition().x(), 0., 1e-9 * mm);
    }},
    {"event_loop_fatal_mt", [] {   // must abort: registered with bbrsim_add_abort_test
      WriteBadRowDataset("BadRowCrack");
      G4RunManager* rm = BuildEventLoop(G4RunManagerType::TaskingOnly, "BadRowCrack");
      rm->Initialize();
      rm->BeamOn(1);
      // Not reached when the worker aborts. If BeamOn returns, the case exits 0, and
      // ExpectAbort.cmake fails it: the test requires the abort itself.
      std::printf("BeamOn returned: the worker's BBR013 did not abort the process\n");
    }},
    {"missing_sidecar_stops_before_events", [] {
      hfssfix::WriteDataset(DataRoot(), "Bare_500GHz", hfssfix::Mini500(), false);
      G4RunManager* rm = BuildEventLoop(G4RunManagerType::SerialOnly, "Bare");
      rm->Initialize();
      ExpectG4Exception("BBR024", [&] { rm->BeamOn(1); }, "BBRDatasetSidecar",
                        "Bare_500GHz.dataset.json: cannot be opened");
      CHECK(G4StateManager::GetStateManager()->GetCurrentState() != G4State_EventProc);
    }},
    {"section_wider_than_crack_stops_before_events", [] {
      auto p = hfssfix::SidecarFrom("Wide_500GHz", hfssfix::Mini500());
      p.zHalf = 3e-5;   // 30 um, beyond the crack's 26 um half-gap
      hfssfix::WriteDataset(DataRoot(), "Wide_500GHz", hfssfix::Mini500(), false);
      bbrtest::WriteFile(DataRoot() / "waveguides" / "Wide_500GHz.dataset.json", hfssfix::SidecarJson(p));
      G4RunManager* rm = BuildEventLoop(G4RunManagerType::SerialOnly, "Wide");
      rm->Initialize();
      ExpectG4Exception("BBR025", [&] { rm->BeamOn(1); }, "BBRDatasetSidecar",
                        "does not fit strictly inside crack volume Wide");
      CHECK(G4StateManager::GetStateManager()->GetCurrentState() != G4State_EventProc);
    }},
  });
}
