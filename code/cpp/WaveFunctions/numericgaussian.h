#pragma once
#include "simplegaussian.h"

class NumericGaussian : public SimpleGaussian {
  /* A simple extension of SimpleGaussian that calculates the Laplacian
     of the wavefunction numerically instead of using an analytic expression
  */

public:
  NumericGaussian(class System* system, double alpha) : SimpleGaussian(system, alpha) { }
  double computeDoubleDerivative(std::vector<class Particle*> particles);

};
