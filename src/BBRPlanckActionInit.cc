#include "BBRPlanckActionInit.hh"
#include "BBRThermalPGA.hh"

void BBRPlanckActionInit::Build() const
{
  SetUserAction(new BBRThermalPGA());
}
