#include "numericgaussian.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
using namespace std;

double NumericGaussian::computeDoubleDerivative(std::vector<class Particle*> particles){

  int numParticles  = m_system->getNumberOfParticles();
  int numDimensions = m_system->getNumberOfDimensions();

  double h       = 1e-5;
  double minus2h = -2*h;

  // Calculate the double derivative, see Drafts folder on github for derivation of formula
  double doubleDerivative = 0, term = 0;
  for (int i = 0; i < numParticles; i++){

    for (int j = 0; j < numDimensions; j++){
      particles[i]->adjustPosition(h, j);
      term +=evaluate(particles, i);
      particles[i]->adjustPosition(minus2h, j);
      term +=evaluate(particles, i);
      particles[i]->adjustPosition(h, j);
    }
    term/=evaluate(particles,i);
    doubleDerivative += term;
    term = 0;
  }
  doubleDerivative -= 6*numParticles;
  doubleDerivative /= h*h;
  return doubleDerivative;
}
