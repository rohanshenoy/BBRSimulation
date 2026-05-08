/// \file BBRSim.cc
/// \brief Main program for BBRSimulation
///        Adapted from Geant4 optical/OpNovice2 example.

#include "ActionInitialization.hh"
#include "BBRCrackDetectorConstruction.hh"
#include "BBRReflectanceDetectorConstruction.hh"
#include "BBRDiffractionActionInit.hh"
#include "BBSimPhysics.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "SteppingVerbose.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if(argc == 1) ui = new G4UIExecutive(argc, argv);

  auto steppingVerbose = new SteppingVerbose;

  auto runManager = G4RunManagerFactory::CreateRunManager();

  // Use the crack diffraction geometry for diffraction*.mac and validation.mac.
  G4String macName      = argc > 1 ? G4String(argv[1]) : G4String("");
  bool isBBR            = macName.find("diffraction") != G4String::npos
                       || macName.find("validation")  != G4String::npos;
  bool isReflectance    = macName.find("reflectance") != G4String::npos;

  G4VUserDetectorConstruction* detector;
  if      (isBBR)         detector = new BBRCrackDetectorConstruction();
  else if (isReflectance) detector = new BBRReflectanceDetectorConstruction();
  else                    detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  auto opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);

  // BBSimPhysics must be registered after G4OpticalPhysics so that
  // G4OpBoundaryProcess already exists in the process manager when
  // BBSimPhysics::ConstructProcess() runs.
  physicsList->RegisterPhysics(new BBSimPhysics());

  runManager->SetUserInitialization(physicsList);

  if (isBBR) {
    G4double gunZ = macName.find("crack2") != G4String::npos ? 3.*mm : 0.;
    runManager->SetUserInitialization(new BBRDiffractionActionInit(gunZ));
  } else if (isReflectance) {
    // Gun default: (−20 mm, 0, 0) along +x — fires at Cu face at x=0.
    runManager->SetUserInitialization(new BBRDiffractionActionInit(0.));
  } else {
    runManager->SetUserInitialization(new ActionInitialization());
  }

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(ui)
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;
  delete steppingVerbose;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
