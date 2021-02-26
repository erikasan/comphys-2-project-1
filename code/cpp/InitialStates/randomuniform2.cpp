#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"
#include "randomuniform2.h"

using std::cout;
using std::endl;

RandomUniformMinDist::RandomUniformMinDist(System* system, int numberOfDimensions, int numberOfParticles, double minDist) :
        RandomUniform(system,numberOfDimensions,numberOfParticles) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    m_minDist=minDist;
    setupInitialState();
}

void RandomUniformMinDist::setupInitialState() {
    double randnumb;
    int counter=0; //How often a random placement has gone wrong
    double multiplier=1; //The number to multiply the random placing with (so we get a larger number)
    bool accepted; //If the particle can be placed
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();
        accepted=true;
        replace: //Label for replacing the particle
        for (int j=0; j < m_numberOfDimensions; j++) {
              //Create a random double between 0 and 1*multiplier to place the particle, and places the particle there
              randnumb=m_system->getRandomEngine()->nextDouble()*multiplier;
              position.push_back(randnumb);
          }
          for(int k=0; k<i;k++){ //For each particle, check wether the distance to the other particles is sufficiently large
            if(distance(m_particles[k]->getPosition(),position)<=m_minDist){
              accepted=false;
              break;
            }
        }
        if(accepted){ //If the move is accepted (if the distance is sufficiently large to all particles)
          m_particles.push_back(new Particle());
          m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
          m_particles.at(i)->setPosition(position);
        }
        else{
          /*
          Clear the positions, increase the number of wrong placements by one,
          if the misplacements ar consistent, increase the multiplier,
          go back to placing the particle
          */
          position.clear();
          counter++;
          if (counter>=10){
            multiplier*=1.1;
            counter=0;
          }
          accepted=true;
          goto replace;
        }
    }
}
double RandomUniformMinDist::distance(std::vector<double> part1, std::vector<double> part2){
  double val=0;
  for (int i=0; i<m_numberOfDimensions;i++){
    val+=(part1[i]-part2[i])*(part1[i]-part2[i]);
  }
  return sqrt(val);
}
