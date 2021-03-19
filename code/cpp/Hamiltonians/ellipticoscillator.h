#pragma once
#include "hamiltonian.h"
#include "harmonicoscillator.h"
#include <vector>

class EllipticOscillator : public HarmonicOscillator {
public:
  //using HarmonicOscillator:HarmonicOscillator; //same constructor type
  EllipticOscillator(System* system, double omega) : HarmonicOscillator(system, omega) {};
  double computeLocalPotentialEnergy(std::vector<Particle*> particles) override;
  // The local potential energy is different, the kinetic energy is still only the derivative
private:
  double gamma=2.82843;
};
