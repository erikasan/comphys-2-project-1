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

    //For each particle
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();
        accepted=true;
        replace: //Label for retrying

        //The for loop suggests a new position
        for (int j=0; j < m_numberOfDimensions; j++) {
              //Create a random double between 0 and 1*multiplier to place the particle, and places the particle there
              randnumb=(m_system->getRandomEngine()->nextDouble()-0.5)*multiplier;
              position.push_back(randnumb);
          }
        //It is checked if the suggested position has the relevant distance from the other particles to be legally placed
        for(int k=0; k<i;k++){
          if(distance(m_particles[k]->getPosition(),position)<=m_minDist){
            accepted=false;
            break;
          }
        }
        //If the move is accepted (if the distance is sufficiently large to all particles), create an new particle
        if(accepted){
          m_particles.push_back(new Particle());
          m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
          m_particles.at(i)->setPosition(position);
        }

        /*
        Clear the positions, increase the number of wrong placements by one,
        if the misplacements are consistent, increase the multiplier (so that new suggestions are further away)
        go back to placing the particle
        */
        else{

          position.clear();
          counter++;
          if (counter>=10){ //If more than 10 placing suggestions have failed
            multiplier*=1.1;
            counter=0; //10 more to fail to be multiplied by 1.1 again
          }
          accepted=true; //Default value is true
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
