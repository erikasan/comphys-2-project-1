#pragma once
#include "randomuniform.h"
class RandomUniformMinDist : public RandomUniform {
public:
    RandomUniformMinDist(System* system, int numberOfDimensions, int numberOfParticles, double minDist);
    void setupInitialState();
private:
  double m_minDist; //The minimal distance upon creation
  double distance(std::vector<double> part1, std::vector<double> part2); //Returns the absolute distance between two particles
};
