#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "../sampler.h"
#include <iostream>
#include <iomanip>
#include "complexfunction.h"

using namespace std;

ComplexFunction::ComplexFunction(System* system, double alpha, double beta_param, double a_param) : WaveFunction(system) {
    assert(alpha >= 0);
    beta = beta_param;
    a = a_param;
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

void ComplexFunction::printMatrix(double **A, int n){ //Simple print matrix function, needed for debugging
  cout << std::setprecision(3) << fixed;
  for (int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
}
double ** ComplexFunction::createNNMatrix(int n){ //Create matrix function, needed for storing distances
  double** A;
  A = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[n];
  }
  for (int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      A[i][j] = 0.0;
    }
  }
  return A;
}

void ComplexFunction::initiateDistances(std::vector<class Particle*> particles){
  //Initate distance matrix for particles, as these are used a lot
  //The matrix is symmetric, but the algorithms get easier by double-storing distances
  particle_distances_absolute=createNNMatrix(m_system->getNumberOfParticles());
  for (int i = 0; i < m_system->getNumberOfParticles(); i++){
      for (int j = 0; j < m_system->getNumberOfParticles(); j++){
        particle_distances_absolute[i][j] = distance(particles[i]->getPosition(), particles[j]->getPosition());
      }
  }

}

void ComplexFunction::updateDistances(std::vector<class Particle*> particles,int particle_id){
  //Update distances in the matrix  from all particles to the particle particle_id
  double distanceval = 0;
  std::vector<double> particlePosition=particles[particle_id]->getPosition();
  for (int i = 0; i < m_system->getNumberOfParticles(); i++){
    distanceval = distance(particlePosition,particles[i]->getPosition());
    particle_distances_absolute[i][particle_id] = distanceval;
    particle_distances_absolute[particle_id][i] = distanceval;
  }
}

double ComplexFunction::distance(std::vector<double> part1, std::vector<double> part2){
  //Distance between two vectors
  double val = 0;
  for (int i = 0; i < m_system->getNumberOfDimensions(); i++){
    val += (part1[i] - part2[i])*(part1[i] - part2[i]);
  }
  return sqrt(val);
}
double ComplexFunction::u(double r){
  return log(1 - a/r);
}
double ComplexFunction::uder(double r){
  return -a/(a*r - r*r);
}
double ComplexFunction::uderder(double r){
  return a*(a - 2*r)/(r*r*(a - r)*(a - r));
}

double ComplexFunction::evaluate(std::vector<class Particle*> particles) {
    //Evaluate the wave function for a fixed set of particle positions
     double alpha = m_parameters[0];
     double total_radius = 0;
     int dimension = m_system->getNumberOfDimensions();
     int numberOfParticles = m_system->getNumberOfParticles();
     for (int i = 0; i < numberOfParticles; i++){
       std::vector<double> position = particles[i]->getPosition();
       for (int j=0; j < dimension-1; j++){
         total_radius += position[j]*position[j];
       }
       total_radius += position[dimension-1]*position[dimension-1]*beta; //Add beta part for the z coordinate
     }
     double prod_g = exp(-alpha*total_radius); // The g-part of the wave function
     double prod_f = 1;
     for (int i = 0; i < numberOfParticles; i++){ // For each particle
       for (int j = i+1; j < numberOfParticles; j++){ // For all particles with a larger value
         if(particle_distances_absolute[i][j] < a){
           return 0; //If the particle distance is illegal (less than a)
         }
         prod_f*=(1 - a/particle_distances_absolute[i][j]);
       }
     }
     return prod_g*prod_f;
}

double ComplexFunction::evaluate(std::vector<class Particle*> particles, int particle_id) {
    //As the Wave function is partially seperable, evaluates only the part belonging to particle "particle_id"
    //Partially seperable in the sense that there is a part of the wave function that is independent when only moving 1 particle
    double alpha = m_parameters[0];
    double total_radius = 0;
    int dimension = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();
    std::vector<double> position = particles[particle_id]->getPosition();
    double prod_f = 1;


    //Calculate f-part of wave function
    for (int i = 0; i < particle_id; i++){ //For all particles with a lower id
      if(particle_distances_absolute[particle_id][i] < a){
        return 0;
      }
      prod_f *= (1 - a/particle_distances_absolute[particle_id][i]);
    }
    for (int i = particle_id+1; i < numberOfParticles; i++){ //For all particles with a larger id
      if(particle_distances_absolute[particle_id][i] < a){
        return 0;
      }
      prod_f *= (1 - a/particle_distances_absolute[particle_id][i]);
    }


    //Calculate G-part of Wave function
    for (int j=0; j < dimension-1; j++){
      total_radius += position[j]*position[j];
    }
    total_radius += position[dimension-1]*position[dimension-1]*beta;
    double prod_g = exp(-alpha*total_radius);

    return prod_g*prod_f;
}

double ComplexFunction::computeDoubleDerivative(std::vector<class Particle*> particles) {
      //Compute the local energy. Sadly, this is not a seperable WaveFunction
      // So it is really expensive to calculate.
      double total_energy = 0;
      double total_radius_withbeta = 0; //This is x^2+y^2+beta^2z^2
      double alpha = m_parameters[0];
      int numberOfParticles = m_system->getNumberOfParticles();
      int num_dim = m_system->getNumberOfDimensions();
      double third_sum_temp = 0;
      double prefactor;
      std::vector<double> distance_vec = std::vector<double>(num_dim,0.0); //Reusable vector, distance between particles
      std::vector<double> temp = std::vector<double>(num_dim,0.0); //Temporary vector for calculations
      for (int k = 0; k < numberOfParticles; k++){
        std::vector<double> position = particles[k]->getPosition();
        std::vector<double> first_sum_vector = position;
        first_sum_vector[num_dim-1] *= beta; //multiply z-parameter by beta

        //Calculate the energy of the non-interacting wave function
        for (int j = 0; j < num_dim-1; j++){
          total_radius_withbeta += position[j]*position[j];
        }
        total_radius_withbeta += position[num_dim-1]*position[num_dim-1]*beta*beta;


        for (int i = 0; i < k; i++){
          //The prefactor in the first sum which is reused for each dimension
          prefactor = -a/((a - particle_distances_absolute[k][i])*particle_distances_absolute[k][i]*particle_distances_absolute[k][i]);
          for (int l = 0; l < num_dim; l++){
            distance_vec[l] = position[l]-particles[i]->getPosition()[l]; //Distance between particle [i] and particle [k] as vector
            temp[l] += prefactor*distance_vec[l];
          }
          third_sum_temp = a/(particle_distances_absolute[k][i]*(a - particle_distances_absolute[k][i])); //The value for the third sum
          total_energy += -(third_sum_temp*third_sum_temp);
        }

        //Repetition of the previous loop for the particles with higher index than k (so k itself isn't included)
        for (int i = k+1; i < numberOfParticles; i++){
          prefactor = -a/((a - particle_distances_absolute[k][i])*particle_distances_absolute[k][i]*particle_distances_absolute[k][i]);
          for (int l = 0; l < num_dim; l++){
            distance_vec[l] = position[l]-particles[i]->getPosition()[l]; //Distance between particle [i] and particle [k] as vector
            temp[l] += prefactor*distance_vec[l];
          }
          third_sum_temp = a/(particle_distances_absolute[k][i]*(a-particle_distances_absolute[k][i])); //The value for the third sum
          total_energy += -(third_sum_temp*third_sum_temp);
        }


        for(int l = 0; l < num_dim; l++){
          total_energy += temp[l]*temp[l];
          total_energy += -4*alpha*temp[l]*first_sum_vector[l];
        }
        std::fill(temp.begin(), temp.end(), 0.0);

      }
      total_energy += 4*alpha*alpha*total_radius_withbeta;
      total_energy -= numberOfParticles*((2*num_dim - 2)*alpha + 2*alpha*beta);
      return total_energy;
}

std::vector<double> ComplexFunction::quantumForce(std::vector<class Particle*> particles, int particle_id){
     double alpha = m_parameters[0];
     int num_dim = m_system->getNumberOfDimensions();
     std::vector<double> qForce = std::vector<double>();
     qForce.reserve(num_dim);
     std::vector<double> temp = std::vector<double>(num_dim, 0.0);
     std::vector<double> particle_position = particles[particle_id]->getPosition();
     for(int l = 0; l < num_dim; l++){
       temp[l] += particle_position[l]*(-4*alpha);
     }
     temp[num_dim-1] *= beta;

     double prefactor;
     double dist;
     std::vector<double> distance_vec = std::vector<double>(num_dim,0.0);
     for (int i = 0; i < particle_id; i++){
       dist = particle_distances_absolute[i][particle_id];
       prefactor = 2*a/((dist*dist)*(a - dist));
       for(int l = 0; l < num_dim; l++){
         distance_vec[l] = particle_position[l]-particles[i]->getPosition()[l];
         temp[l] += -prefactor*distance_vec[l];
       }
     }
     for (int i = particle_id+1; i < m_system->getNumberOfParticles(); i++){
       dist = particle_distances_absolute[i][particle_id];

       prefactor =2*a/((dist*dist)*(a - dist));
       for(int l = 0; l < num_dim; l++){
         distance_vec[l] = particle_position[l]-particles[i]->getPosition()[l];
         temp[l] += -prefactor*distance_vec[l];
       }
     }
     for (int i = 0; i < num_dim; i++){
       qForce.push_back(temp[i]);
     }
     return qForce;
}

std::vector<double> ComplexFunction::quantumForce(std::vector<class Particle*> particles){
  //Unused, but might be useful later.
  int num_dim = m_system->getNumberOfDimensions();
  int num_part = m_system->getNumberOfParticles();
  std::vector<double> qForce = std::vector<double>();
  std::vector<double> qForce_oneparticle = std::vector<double>();
  qForce.reserve(num_dim*num_part);
  qForce_oneparticle.reserve(num_dim);
  for (int i = 0; i < num_part; i++){
    qForce_oneparticle = quantumForce(particles,i);
    for (int j = 0; j < num_dim; j++){
      qForce[num_dim*i+j] = qForce_oneparticle[j];
    }
  }
  return qForce;
}

double ComplexFunction::totalRadius(std::vector<class Particle*> particles){
  double total_radius = 0;
  std::vector<double> position;
  position.reserve(3);

  int num_part = m_system->getNumberOfParticles();
  for (int i = 0; i < num_part; i++){
    position = particles[i]->getPosition();
    total_radius += position[0]*position[0] + position[1]*position[1] + beta*position[2]*position[2];
  }

  return total_radius;
}

double ComplexFunction::localEnergyTotalRadius(std::vector<class Particle*> particles, double localEnergy){
  double total_radius = ComplexFunction::totalRadius(particles);
  return localEnergy*total_radius;
}

void ComplexFunction::sample(std::vector<class Particle*> particles, double localEnergy){
  m_av_total_radius += ComplexFunction::totalRadius(particles);
  m_av_local_energy_total_radius += ComplexFunction::localEnergyTotalRadius(particles, localEnergy);
  return;
}

void ComplexFunction::computeAverages(double steps){
  m_av_total_radius /= steps;
  m_av_local_energy_total_radius /= steps;
  return;
}

void ComplexFunction::gradientDescent(){
  double alphaOld, alphaNew, energy, localEnergyGradient;

  alphaOld            = m_parameters[0];
  energy              = m_system->getSampler()->getEnergy();
  localEnergyGradient = 2*(energy*m_av_total_radius - m_av_local_energy_total_radius);

  alphaNew = alphaOld - m_learningRate*localEnergyGradient;

  if (abs(alphaNew - alphaOld) < m_tol){
    m_system->stopGradientDescent();
  }

  m_parameters[0] = alphaNew;

  m_av_total_radius = m_av_local_energy_total_radius = 0;
  return;
}
