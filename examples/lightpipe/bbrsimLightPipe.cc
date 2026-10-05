/// \file bbrsimLightPipe.cc
/// \brief Entry point for the light-pipe geometry (parametric or CAD).

#include "LightPipeDetectorConstruction.hh"
#include "LightPipeActionInitialization.hh"
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

  auto* runManager = G4RunManagerFactory::CreateRunManager();
  BBRConfigManager::Instance();

  runManager->SetUserInitialization(new LightPipeDetectorConstruction());

  auto* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->RegisterPhysics(new BBSimPhysics());
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new LightPipeActionInitialization());

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
