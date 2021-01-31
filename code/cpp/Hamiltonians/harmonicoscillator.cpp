#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}
double HarmonicOscillator::computeLocalKineticEnergy(std::vector<Particle*> particles){
  double kineticEnergy   = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
  return kineticEnergy;
}
double HarmonicOscillator::computeLocalPotentialEnergy(std::vector<Particle*> particles){
  double potentialEnergy =0;
  for (int i=0; i< m_system->getNumberOfParticles();i++){
    potentialEnergy+=particles[i]->getRadiussquared();
  }
  potentialEnergy*=0.5*m_omega*m_omega;
  return potentialEnergy;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */
    return computeLocalPotentialEnergy(particles) + computeLocalKineticEnergy(particles);
}
