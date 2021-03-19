#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "ellipticoscillator.h"
using std::cout;
using std::endl;
double EllipticOscillator::computeLocalPotentialEnergy(std::vector<Particle*> particles){
  double potentialEnergy = 0;
  int dimension=m_system->getNumberOfDimensions();
  for (int i = 0; i < m_system->getNumberOfParticles(); i++){ //for each particle
    std::vector<double> position=particles[i]->getPosition();

    for (int j=0;j<dimension-1;j++){ //For each dimension
      potentialEnergy += position[j]*position[j];
    }
    potentialEnergy+=position[dimension-1]*position[dimension-1]*gamma*gamma; //Gamma factor for z / last dimenion
  }
  potentialEnergy *= 0.5;
  return potentialEnergy;
}
