#ifndef BBRPlanckActionInit_hh
#define BBRPlanckActionInit_hh

#include "G4VUserActionInitialization.hh"

class BBRPlanckActionInit : public G4VUserActionInitialization {
 public:
  BBRPlanckActionInit()  = default;
  ~BBRPlanckActionInit() override = default;
  void Build() const override;
};

#endif
