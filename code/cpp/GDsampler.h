#pragma once
#include "sampler.h"

class GDsampler : public Sampler {

public:
  GDsampler(class System* system) : Sampler(system) { }
  void gdsampler(std::vector<class Particle*> particles, double localEnergy);
};
