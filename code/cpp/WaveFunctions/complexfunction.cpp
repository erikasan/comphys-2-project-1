
#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>
#include "complexfunction.h"
using namespace std;
ComplexFunction::ComplexFunction(System* system, double alpha) : WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}
double ** ComplexFunction::createNNMatrix(int n){
  double** A;
  A = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[n];
  }
  for (int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A[i][j]=0.0;
    }
  }
  return A;
}
void ComplexFunction::calculateParticleDistances(std::vector<class Particle*> particles){
  for (int i=0; i<m_system->getNumberOfParticles();i++){
      for (int j=0; j<m_system->getNumberOfParticles();j++){
        particle_distances_absolute[i][j]=distance(particles[i]->getPosition(),particles[j]->getPosition());
      }
  }
}
void ComplexFunction::updateParticleDistances(std::vector<class Particle*> particles,int particle_id){
  double distanceval=0;
  std::vector<double> particlePosition=particles[particle_id]->getPosition();
  for (int i=0; i<m_system->getNumberOfParticles();i++){
    distanceval=distance(particlePosition,particles[i]->getPosition());
    particle_distances_absolute[i][particle_id]=distanceval;
    particle_distances_absolute[particle_id][i]=distanceval;
  }
}
double ComplexFunction::vectorProduct(std::vector<double> part1, std::vector<double> part2){
  double val=0;
  for (int i=0; i<m_system->getNumberOfDimensions();i++){
    val+=part1[i]*part2[i];
  }
  return val;
}
double ComplexFunction::distance(std::vector<double> part1, std::vector<double> part2){
  double val=0;
  for (int i=0; i<m_system->getNumberOfDimensions();i++){
    val+=part1[i]*part1[i]-part2[i]*part2[i];
  }
  return sqrt(val);
}
std::vector<double> ComplexFunction:: distance_vector(std::vector<double> part1, std::vector<double> part2){
  int num_dim=m_system->getNumberOfDimensions();
  std::vector<double> distance_vec = std::vector<double>();
  distance_vec.reserve(num_dim);
  for (int i=0; i<m_system->getNumberOfDimensions();i++){
    distance_vec.push_back(part1[i]-part2[i]);
  }
  return distance_vec;
}
double ComplexFunction::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
     double alpha = m_parameters[0];
     double total_radius = 0;
     int dimension=m_system->getNumberOfDimensions();
     for (int i = 0; i < m_system->getNumberOfParticles(); i++){
       std::vector<double> position=particles[i]->getPosition();
       for (int j=0;j<dimension-1;j++){
         total_radius += position[j]*position[j];
       }
       total_radius+=position[dimension-1]*position[dimension-1]*beta;
     }
     for(Particle *particle : particles){
       total_radius += particle->getRadiussquared();
     }
     double prod_g=exp(-alpha*total_radius);
     double prod_f=1;
     for (int i=0; i<m_system->getNumberOfParticles(); i++){
       for (int j=i+1; j<m_system->getNumberOfParticles(); j++){
         prod_f*=(1-a/particle_distances_absolute[i][j]);
       }
     }
     return prod_g*prod_f;
}

double ComplexFunction::evaluate(std::vector<class Particle*> particles, int particle_id) {
    //As the Wave function is seperable, evaluates only the part belonging to particle "particle_id"
     double alpha = m_parameters[0];
     double total_radius = 0;
     int dimension=m_system->getNumberOfDimensions();
     std::vector<double> position=particles[particle_id]->getPosition();
     for (int j=0;j<dimension-1;j++){
       total_radius += position[j]*position[j];
     }
     total_radius+=position[dimension-1]*position[dimension-1]*beta;
     double radius = particles[particle_id]->getRadiussquared();
     double g=exp(-alpha*radius);
     double prod_f=1;
     for (int i=0;i<particle_id;i++){
       prod_f*=(1-a/particle_distances_absolute[particle_id][i]);
     }
     for (int i=particle_id+1;i<m_system->getNumberOfParticles();i++){
       prod_f*=(1-a/particle_distances_absolute[particle_id][i]);
     }
     return g*prod_f;
}


std::vector<double> ComplexFunction::quantumForce(std::vector<class Particle*> particles, int particle_id){
     double alpha = m_parameters[0];
     int num_dim=m_system->getNumberOfDimensions();
     std::vector<double> qForce = std::vector<double>();
     qForce.reserve(num_dim);
     std::vector<double> particle_position=particles[particle_id]->getPosition();
     for (int j=0; j < num_dim-1; j++) {
         qForce.push_back(-4*alpha*particle_position[j]);
     }
     qForce.push_back(-4*alpha*beta*particle_position[num_dim-1]);
     double prefactor;
     double dist;
     std::vector<double> distance_vec = std::vector<double>();
     for (int i=0;i<particle_id;i++){
       dist=particle_distances_absolute[particle_id][i];
       prefactor=2*a/((dist*dist)*(a-dist));
       distance_vec=distance_vector(particle_position,particles[i]->getPosition());
       for (int j=0; j < num_dim; j++) {
           qForce[j]-=prefactor*distance_vec[i];
       }
     }
     for (int i=0;i<particle_id;i++){
       dist=particle_distances_absolute[particle_id][i];
       prefactor=2*a/((dist*dist)*(a-dist));
       distance_vec=distance_vector(particle_position,particles[i]->getPosition());
       for (int j=0; j < num_dim; j++) {
           qForce[j]-=prefactor*distance_vec[i];
       }
     }
     return qForce;
}
/*
double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     *
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
*/
