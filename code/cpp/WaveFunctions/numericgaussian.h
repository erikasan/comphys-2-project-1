#pragma once
#include "simplegaussian.h"

class NumericGaussian : public SimpleGaussian {
public:
  double computeDoubleDerivative(std::vector<class Particle*> particles);

};
