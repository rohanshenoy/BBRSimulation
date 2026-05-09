#ifndef BBRTestActionInit_hh
#define BBRTestActionInit_hh
#include "G4VUserActionInitialization.hh"

class BBRTestActionInit : public G4VUserActionInitialization {
public:
  BBRTestActionInit()           = default;
  ~BBRTestActionInit() override = default;
  void Build() const override;
};
#endif
