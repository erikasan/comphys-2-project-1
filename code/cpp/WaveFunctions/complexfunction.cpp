#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>
#include <iomanip>
#include "complexfunction.h"
#include <armadillo>
using namespace std;
using namespace arma;
ComplexFunction::ComplexFunction(System* system, double alpha, double beta_param, double a_param) : WaveFunction(system) {
    assert(alpha >= 0);
    beta=beta_param;
    a=a_param;
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}
void ComplexFunction::printMatrix(double **A, int n){
  cout << std::setprecision(3) << fixed;
  for (int i=0; i<n;i++){
    for(int j=0;j<n;j++){
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
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
void ComplexFunction::initiateDistances(std::vector<class Particle*> particles){
  particle_distances_absolute=createNNMatrix(m_system->getNumberOfParticles());
  for (int i=0; i<m_system->getNumberOfParticles();i++){
      for (int j=0; j<m_system->getNumberOfParticles();j++){
        particle_distances_absolute[i][j]=distance(particles[i]->getPosition(),particles[j]->getPosition());
      }
  }
  printMatrix(particle_distances_absolute,m_system->getNumberOfParticles());
}
void ComplexFunction::updateDistances(std::vector<class Particle*> particles,int particle_id){
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
    val+=(part1[i]-part2[i])*(part1[i]-part2[i]);
  }
  return sqrt(val);
}
double ComplexFunction::u(double r){
  return log(1-a/r);
}
double ComplexFunction::uder(double r){
  return -a/(a*r-r*r);
}
double ComplexFunction::uderder(double r){
  return a*(a-2*r)/(r*r*(a-r)*(a-r));
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
     int numberOfParticles=m_system->getNumberOfParticles();
     for (int i = 0; i < numberOfParticles; i++){
       std::vector<double> position=particles[i]->getPosition();
       for (int j=0;j<dimension-1;j++){
         total_radius += position[j]*position[j];
       }
       total_radius+=position[dimension-1]*position[dimension-1]*beta;
     }
     double prod_g=exp(-alpha*total_radius);
     double prod_f=1;
     for (int i=0; i<numberOfParticles; i++){
       for (int j=i+1; j<numberOfParticles; j++){
         if(particle_distances_absolute[i][j]<a){
           return 0;
         }
         prod_f*=(1-a/particle_distances_absolute[i][j]);
       }
     }
     return prod_g*prod_f;
}

double ComplexFunction::evaluate(std::vector<class Particle*> particles, int particle_id) {
    //As the Wave function is partially seperable, evaluates only the part belonging to particle "particle_id"
    double alpha = m_parameters[0];
    double total_radius = 0;
    int dimension=m_system->getNumberOfDimensions();
    int numberOfParticles=m_system->getNumberOfParticles();
    std::vector<double> position=particles[particle_id]->getPosition();
    double prod_f=1;


    //Calculate f-part of wave function
    for (int i=0;i<particle_id;i++){
      if(particle_distances_absolute[particle_id][i]<a){
        return 0;
      }
      prod_f*=(1-a/particle_distances_absolute[particle_id][i]);
    }
    for (int i=particle_id+1;i<numberOfParticles;i++){
      if(particle_distances_absolute[particle_id][i]<a){
        return 0;
      }
      prod_f*=(1-a/particle_distances_absolute[particle_id][i]);
    }


    //Calculate G-part of Wave function
    for (int j=0;j<dimension-1;j++){
      total_radius += position[j]*position[j];
    }
    total_radius+=position[dimension-1]*position[dimension-1]*beta;
    double prod_g=exp(-alpha*total_radius);

    return prod_g*prod_f;
}

/*
double ComplexFunction::computeDoubleDerivative(std::vector<class Particle*> particles){
  double total_energy=0;
  double total_radius_withbeta = 0; //This is x^2+y^2+beta^2z^2
  double alpha = m_parameters[0];
  int numberOfParticles = m_system->getNumberOfParticles();
  int num_dim  = m_system->getNumberOfDimensions();
  double third_sum_temp=0;
  double prefactor;
  vec distance_vec;
  double rad_alt=0;
  double distance_abs;
  vec  temp(num_dim); //Vector of length num_dim with all 0s
  temp.fill(0.0);
  for (int k = 0; k < numberOfParticles; k++){
    vec position=conv_to<vec>::from(particles[k]->getPosition());
    rad_alt=position[0]*position[0]+position[1]*position[1]+position[2]*position[2]*beta*beta;
    total_energy+=4*alpha*alpha*rad_alt-4*alpha-2*alpha*beta;
    vec first_sum_vector=vec(position);
    first_sum_vector(2)*=beta; //multiply z-parameter by beta
    for (int i=0;i<numberOfParticles;i++){
      if (k==i){
        continue;
      }
      distance_vec=position-conv_to<vec>::from(particles[i]->getPosition());
      distance_abs=particle_distances_absolute[i][k];
      total_energy+=uder(distance_abs)/distance_abs*(-4)*alpha*dot(first_sum_vector,distance_vec);
    }
    for (int i=0;i<numberOfParticles;i++){
        if (k==i){
          continue;
        }
        distance_vec=position-conv_to<vec>::from(particles[i]->getPosition());
        distance_abs=particle_distances_absolute[k][i];
        temp+=distance_vec*(uder(distance_abs)/distance_abs);
    }
    total_energy+=dot(temp,temp);
    temp.fill(0.0);
    for (int i=0;i<numberOfParticles;i++){
      if (k==i){
        continue;
      }
      distance_abs=particle_distances_absolute[i][k];
      total_energy+=uderder(distance_abs)+2/distance_abs*uder(distance_abs);
    }

  }
  return total_energy;
}
*/
double ComplexFunction::computeDoubleDerivative(std::vector<class Particle*> particles) {
      double total_energy=0;
      double total_radius_withbeta = 0; //This is x^2+y^2+beta^2z^2
      double alpha = m_parameters[0];
      int numberOfParticles = m_system->getNumberOfParticles();
      int num_dim  = m_system->getNumberOfDimensions();
      double third_sum_temp=0;
      double prefactor;
      vec distance_vec;
      vec  temp(num_dim); //Vector of length num_dim with all 0s
      temp.fill(0.0);
      for (int k = 0; k < numberOfParticles; k++){
        vec position=conv_to<vec>::from(particles[k]->getPosition());
        vec first_sum_vector=vec(position);
        first_sum_vector(num_dim-1)*=beta; //multiply z-parameter by beta

        //Calculate the energy of the non-interacting wave function
        for (int j=0;j<num_dim-1;j++){
          total_radius_withbeta += position[j]*position[j];
        }
        total_radius_withbeta+=position[num_dim-1]*position[num_dim-1]*beta*beta;


        for (int i=0;i<k;i++){
          distance_vec=position-conv_to<vec>::from(particles[i]->getPosition()); //Distance between particle [i] and particle [k] as vector
          prefactor=-a/((a-particle_distances_absolute[k][i])*particle_distances_absolute[k][i]*particle_distances_absolute[k][i]); //The prefactor for the 2. sum
          temp+=prefactor*distance_vec; //Add the vector to the temp_vector (which becomes the second sum)
          third_sum_temp=a/(particle_distances_absolute[k][i]*(a-particle_distances_absolute[k][i])); //The value for the third sum
          total_energy+=-(third_sum_temp*third_sum_temp);
        }

        //Repetition of the previous loop for the particles with higher index than k (so k itself isn't included)
        for (int i=k+1;i<numberOfParticles;i++){
          distance_vec=position-conv_to<vec>::from(particles[i]->getPosition());
          prefactor=-a/((a-particle_distances_absolute[k][i])*particle_distances_absolute[k][i]*particle_distances_absolute[k][i]);
          temp+=prefactor*distance_vec;
          third_sum_temp=a/(particle_distances_absolute[k][i]*(a-particle_distances_absolute[k][i]));
          total_energy+=-(third_sum_temp*third_sum_temp);
        }
        total_energy+=dot(temp,temp);
        total_energy+=-4*alpha*dot(temp,first_sum_vector);
        temp.fill(0.0);
      }
      total_energy+=4*alpha*alpha*total_radius_withbeta;
      total_energy-=numberOfParticles*((2*num_dim-2)*alpha+2*alpha*beta);
      //printf("total_energy %f\n",total_energy);

      return total_energy;
}

std::vector<double> ComplexFunction::quantumForce(std::vector<class Particle*> particles, int particle_id){
     double alpha = m_parameters[0];
     int num_dim=m_system->getNumberOfDimensions();
     std::vector<double> qForce = std::vector<double>();
     qForce.reserve(num_dim);
     vec  temp(num_dim); //Vector of length num_dim with all 0s
     temp.fill(0.0);
     vec particle_position=conv_to<vec>::from(particles[particle_id]->getPosition());
     temp+=particle_position*(-4*alpha);
     temp(num_dim-1)*=beta;
     double prefactor;
     double dist;
     vec distance_vec;
     for (int i=0;i<particle_id;i++){
       dist=particle_distances_absolute[i][particle_id];
       prefactor=2*a/((dist*dist)*(a-dist));
       distance_vec=particle_position-conv_to<vec>::from(particles[i]->getPosition());
       temp+=-prefactor*distance_vec[i];
     }
     for (int i=particle_id+1;i< m_system->getNumberOfDimensions();i++){
       dist=particle_distances_absolute[particle_id][i];
       prefactor=1/dist*uder(dist);
       distance_vec=particle_position-conv_to<vec>::from(particles[i]->getPosition());
       temp+=-prefactor*distance_vec[i];
     }
     for (int i=0;i<num_dim;i++){
       qForce.push_back(temp(i));
     }
     return qForce;
}
std::vector<double> ComplexFunction::quantumForce(std::vector<class Particle*> particles){
  //UNFINISHED!!!!!!!!!!!!
  double alpha = m_parameters[0];
  int num_dim=m_system->getNumberOfDimensions();
  std::vector<double> qForce = std::vector<double>();
  qForce.reserve(num_dim);
  return qForce;
}
