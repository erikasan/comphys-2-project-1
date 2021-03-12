#include "system.h"


class MetropolisLangevin : public System {

public:
  MetropolisLangevin() : System() { }
  MetropolisLangevin(int seed) : System(seed) { }
  bool metropolisStep();
  void runMetropolisSteps(int numberOfMetropolisSteps, bool desire_output);

};
