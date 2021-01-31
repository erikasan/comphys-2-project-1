#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double alpha) : WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
     double alpha = m_parameters[0];
     double total_radius = 0;

     for(Particle *particle : particles){
       total_radius += particle->getRadiussquared();
     }
     return exp(-alpha*total_radius);
}
double SimpleGaussian::evaluate(std::vector<class Particle*> particles, int particle_id) {
    //As the Wave function is seperable, evaluates only the part belonging to particle "particle_id"
     double alpha = m_parameters[0];
     double radius = particles[particle_id]->getRadiussquared();
     return exp(-alpha*radius);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
     double dobdev = 0; // Not needed?
     double total_radius = 0;
     double alpha = m_parameters[0];
     for(Particle *particle : particles){
       total_radius += particle->getRadiussquared();
     }
     int num_part = m_system->getNumberOfParticles();
     int num_dim  = m_system->getNumberOfDimensions();
     return (-2*num_part*num_dim*alpha + 4*alpha*alpha*total_radius);
}

double SimpleGaussian::computeDoubleDerivativeNumerically(std::vector<class Particle*> particles) {
  /* NOTE: This works for a general wavefunction with the only requirement
     that it is separable with respect to the particles (a property utilized to
     make the computation more efficient). The function could therefore
     arguably be moved to a superclass, say the WaveFunction superclass or
     the Hamiltonian superclass.
  */

  int numParticles  = m_system->getNumberOfParticles();
  int numDimensions = m_system->getNumberOfDimensions();

  double h       = 1e-4; // Maybe move this somewhere else?
  double minus2h = -2*h;

  std::vector<class Particle>  tmpParticles;
  std::vector<class Particle*> tmpParticlePointers;
  /* Computing the derivative numerically involves evaluating the WaveFunction
     a stepsize away from the current position. However, as this is not possible
     with the current setup without adjusting the positions of the particles
     using Particle::adjustPosition and then using SimpleGaussian::evaluate,
     a temporary vector of Particle objects (NOT pointers) is required to make
     absolutely certain that the ACTUAL positions of the particles is not changed
     by calculating the derivative.

     SimpleGaussian::evaluate requires a vector of Particle pointers,
     hence a temporary vector of Particle pointers is also declared.
  */

  // Fill tmpParticles with Particle instances and tmpParticlePointers with pointers to Particle
  int i = 0;
  for (Particle *particle : particles){
    tmpParticles.push_back(*particle);
    tmpParticlePointers.push_back(&(tmpParticles[i]));
    i++;
  }

  // Calculate the double derivative, see Drafts folder on github for derivation of formula
  double doubleDerivative = 0;
  for (int i = 0; i < numParticles; i++){
    for (int j = 0; j < numDimensions; j++){
      tmpParticlePointers[i]->adjustPosition(h, j);
      doubleDerivative += SimpleGaussian::evaluate(tmpParticlePointers, i);
      tmpParticlePointers[i]->adjustPosition(minus2h, j);
      doubleDerivative += SimpleGaussian::evaluate(tmpParticlePointers, i);
      tmpParticlePointers[i]->adjustPosition(h, j);
    }
    doubleDerivative /= SimpleGaussian::evaluate(tmpParticlePointers, i);
    doubleDerivative -= 6;
  }
  doubleDerivative /= h*h;
  return doubleDerivative;
}
