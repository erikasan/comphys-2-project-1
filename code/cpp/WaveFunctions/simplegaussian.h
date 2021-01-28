#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction { //Gaussian thus extends WaveFunction
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double evaluate(std::vector<class Particle*> particles, int particle_id);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles, int particle_id);
};
