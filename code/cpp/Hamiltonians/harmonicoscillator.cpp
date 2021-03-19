#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) : Hamiltonian(system) {
  assert(omega > 0);
  m_omega = omega;
}

double HarmonicOscillator::computeLocalKineticEnergy(std::vector<Particle*> particles) {
  double kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
  return kineticEnergy;
}


double HarmonicOscillator::computeLocalPotentialEnergy(std::vector<Particle*> particles) {
  double potentialEnergy = 0;
  for (int i = 0; i < m_system->getNumberOfParticles(); i++){
    potentialEnergy += particles[i]->getRadiussquared();
  }
  potentialEnergy *= 0.5;
  return potentialEnergy;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    //The local energy is just the sum of the potential local energy and the kinetic local energy
    return computeLocalPotentialEnergy(particles) + computeLocalKineticEnergy(particles);
}
