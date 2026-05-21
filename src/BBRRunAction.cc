#include "BBRRunAction.hh"
#include "G4Run.hh"

BBRRunAction::BBRRunAction() : G4UserRunAction() {}

BBRRunAction::~BBRRunAction()
{
  if (fOut.is_open()) fOut.close();
}

void BBRRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "=== BBR Run " << run->GetRunID() << " begin ===" << G4endl;
  if (fOut.is_open()) fOut.close();
  fOut.open("test_output.csv");
  fOut << "event_id,x_mm,y_mm,z_mm,energy_eV,"
          "px_pre,py_pre,pz_pre,px_post,py_post,pz_post,"
          "theta_in_deg,phi_in_deg,"
          "vol_pre,mat_pre,vol_post,mat_post,status,n_reflect\n";
}

void BBRRunAction::EndOfRunAction(const G4Run* run)
{
  if (fOut.is_open()) fOut.close();
  G4cout << "=== BBR Run " << run->GetRunID()
         << " end: test_output.csv written ===" << G4endl;
}
