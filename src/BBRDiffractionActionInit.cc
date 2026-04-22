#include "BBRDiffractionActionInit.hh"

#include "BBRDiffractionPGA.hh"
#include "RunAction.hh"

void BBRDiffractionActionInit::BuildForMaster() const
{
  SetUserAction(new RunAction(nullptr));
}

void BBRDiffractionActionInit::Build() const
{
  SetUserAction(new BBRDiffractionPGA());
  SetUserAction(new RunAction(nullptr));
}
