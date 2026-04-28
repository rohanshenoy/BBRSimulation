#include "BBRDiffractionActionInit.hh"

#include "BBRDiffractionPGA.hh"
#include "RunAction.hh"

BBRDiffractionActionInit::BBRDiffractionActionInit(G4double gunZ_mm)
  : fGunZ_mm(gunZ_mm) {}

void BBRDiffractionActionInit::BuildForMaster() const
{
  SetUserAction(new RunAction(nullptr));
}

void BBRDiffractionActionInit::Build() const
{
  SetUserAction(new BBRDiffractionPGA(fGunZ_mm));
  SetUserAction(new RunAction(nullptr));
}
