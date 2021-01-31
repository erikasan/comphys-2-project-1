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
     int num_dim = m_system->getNumberOfDimensions();
     return (-2*num_part*num_dim*alpha + 4*alpha*alpha*total_radius);
}
