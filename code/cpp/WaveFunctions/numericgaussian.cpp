#include "numericgaussian.h"
#include "../system.h"
#include "../particle.h"



double NumericGaussian::computeDoubleDerivative(std::vector<class Particle*> particles){


  int numParticles  = m_system->getNumberOfParticles();
  int numDimensions = m_system->getNumberOfDimensions();

  double h       = 1e-5; // Maybe move this somewhere else?
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
  double doubleDerivative = 0, term = 0;
  for (int i = 0; i < numParticles; i++){
    for (int j = 0; j < numDimensions; j++){
      particles[i]->adjustPosition(h, j);
      term += SimpleGaussian::evaluate(particles, i);
      particles[i]->adjustPosition(minus2h, j);
      term += SimpleGaussian::evaluate(particles, i);
      //cout << "breaker2 "<< tmpParticlePointers[i]->getPosition().size() << endl; error?
      particles[i]->adjustPosition(h, j);
      /*
      tmpParticlePointers[i]->adjustPosition(h, j);
      term += SimpleGaussian::evaluate(tmpParticlePointers, i);
      tmpParticlePointers[i]->adjustPosition(minus2h, j);
      cout << "breaker1 "<< j << endl;
      term += SimpleGaussian::evaluate(tmpParticlePointers, i);
      //cout << "breaker2 "<< tmpParticlePointers[i]->getPosition().size() << endl; error?
      tmpParticlePointers[i]->adjustPosition(h, j);
    }
    term /= SimpleGaussian::evaluate(tmpParticlePointers, i);
      */
    }
    term /= SimpleGaussian::evaluate(particles, i);
    term -= 6;
    doubleDerivative += term;
    term = 0;
  }
  doubleDerivative /= h*h;
  return doubleDerivative;
}
