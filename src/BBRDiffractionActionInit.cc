#include "BBRDiffractionActionInit.hh"

#include "BBRDiffractionPGA.hh"

BBRDiffractionActionInit::BBRDiffractionActionInit(G4double gunZ_mm)
  : fGunZ_mm(gunZ_mm) {}

void BBRDiffractionActionInit::BuildForMaster() const {}

void BBRDiffractionActionInit::Build() const
{
  SetUserAction(new BBRDiffractionPGA(fGunZ_mm));
}
