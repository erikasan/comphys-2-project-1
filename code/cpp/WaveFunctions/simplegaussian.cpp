#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>
using namespace std;
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
std::vector<double> SimpleGaussian::quantumForce(std::vector<class Particle*> particles){
  std::vector<double> qForce = std::vector<double>();
  double alpha = m_parameters[0];
  qForce.reserve(m_system->getNumberOfDimensions()*m_system->getNumberOfParticles());
  for(int i=0; i< m_system->getNumberOfParticles(); i++){
    std::vector<double> particle_position=particles[i]->getPosition();
    for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
        qForce.push_back(-4*alpha*particle_position[j]);
    }
  }
  return qForce;
}
std::vector<double> SimpleGaussian::quantumForce(std::vector<class Particle*> particles, int particle_id){
     double alpha = m_parameters[0];
     std::vector<double> qForce = std::vector<double>();
     qForce.reserve(m_system->getNumberOfDimensions());
     std::vector<double> particle_position=particles[particle_id]->getPosition();
     for (int j=0; j < m_system->getNumberOfDimensions(); j++) {
         qForce.push_back(-4*alpha*particle_position[j]);
     }
     return qForce;
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
     double total_radius = SimpleGaussian::totalRadius(particles);
     double alpha = m_parameters[0];
     int num_part = m_system->getNumberOfParticles();
     int num_dim  = m_system->getNumberOfDimensions();
     //cout << (-2*num_part*num_dim*alpha + 4*alpha*alpha*total_radius) << endl;
     return (-2*num_part*num_dim*alpha + 4*alpha*alpha*total_radius);
}



double SimpleGaussian::totalRadius(std::vector<class Particle*> particles){
  /* Sum up the squares of the particle distances from the origin.
     This quantity is required for
     - computing the gradient of the average local energy with respect to the variational parameter.
     - computing the Laplacian of the wavefunction.
  */
  double total_radius = 0;
  for (Particle *particle : particles){
    total_radius += particle->getRadiussquared();
  }
  return total_radius;
}

double SimpleGaussian::localEnergyTotalRadius(std::vector<class Particle*> particles, double localEnergy){
  double total_radius = SimpleGaussian::totalRadius(particles);
  return localEnergy*total_radius;
}

void SimpleGaussian::sample(std::vector<class Particle*> particles, double localEnergy){
  m_av_total_radius += SimpleGaussian::totalRadius(particles);
  m_av_local_energy_total_radius += SimpleGaussian::localEnergyTotalRadius(particles, localEnergy);
  return;
}
