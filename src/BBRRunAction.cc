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
  // Open once per session and keep appending: a multi-run macro (e.g. two
  // /run/beamOn with a temperature or gun change in between) must not
  // truncate the previous run's rows. Rows carry run_id to stay separable.
  if (!fOut.is_open()) {
    fOut.open("test_output.csv");
    fOut << "run_id,event_id,x_mm,y_mm,z_mm,energy_eV,"
            "px_pre,py_pre,pz_pre,px_post,py_post,pz_post,"
            "theta_in_deg,phi_in_deg,"
            "vol_pre,mat_pre,vol_post,mat_post,status,n_reflect\n";
  }
}

void BBRRunAction::EndOfRunAction(const G4Run* run)
{
  fOut.flush();
  G4cout << "=== BBR Run " << run->GetRunID()
         << " end: test_output.csv written ===" << G4endl;
}
