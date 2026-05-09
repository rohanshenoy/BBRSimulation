#include "BBRTestActionInit.hh"
#include "BBRTestPGA.hh"
#include "BBRTestSteppingAction.hh"

void BBRTestActionInit::Build() const
{
  SetUserAction(new BBRTestPGA());
  SetUserAction(new BBRTestSteppingAction());
}
