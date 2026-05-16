/// \file BBRSim.cc
/// \brief Main program for BBRSimulation. Adapted from Geant4 OpNovice2.

#include "BBRTestDetectorConstruction.hh"
#include "BBRTestActionInit.hh"
#include "BBRMaterials.hh"
#include "BBSimPhysics.hh"
#include "SteppingVerbose.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include <string>

int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  auto* steppingVerbose = new SteppingVerbose;
  auto* runManager      = G4RunManagerFactory::CreateRunManager();

  // Select Cu wall material from mac filename.
  G4Material* cuMat = nullptr;
  if (argc > 1) {
    std::string mac(argv[1]);
    if (mac.find("of_cu") != std::string::npos)
      cuMat = BBRMaterials::GetOFCopperSerov();
    else if (mac.find("hp_cu") != std::string::npos)
      cuMat = BBRMaterials::GetHPCopperSerov();
    // else nullptr → BBRTestDetectorConstruction defaults to OFHC_Cu
  }

  runManager->SetUserInitialization(new BBRTestDetectorConstruction(cuMat));

  auto* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->RegisterPhysics(new BBSimPhysics());
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new BBRTestActionInit());

  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (ui) {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    UImanager->ApplyCommand(G4String("/control/execute ") + G4String(argv[1]));
  }

  delete visManager;
  delete runManager;
  delete steppingVerbose;
  return 0;
}
