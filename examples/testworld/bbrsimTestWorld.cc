/// \file bbrsimTestWorld.cc
/// \brief Main program for the BBRsim test world (Cu slab with two HFSS crack volumes).

#include "TestWorldDetectorConstruction.hh"
#include "TestWorldActionInitialization.hh"
#include "BBSimPhysics.hh"
#include "BBRConfigManager.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include <cstdlib>

int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  // Multithreaded: G4Analysis merges per-thread ntuples into output/bbr.root.
  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Construct the master config manager + its messenger before any /bbr/
  // command is parsed, so commands work in PreInit and broadcast to workers.
  BBRConfigManager::Instance();

  runManager->SetUserInitialization(new TestWorldDetectorConstruction());

  auto* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->RegisterPhysics(new BBSimPhysics());
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new TestWorldActionInitialization());

  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int commandStatus = 0;
  if (ui) {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    commandStatus = UImanager->ApplyCommand(G4String("/control/execute ") + G4String(argv[1]));
  }

  delete visManager;
  delete runManager;
  return commandStatus == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
