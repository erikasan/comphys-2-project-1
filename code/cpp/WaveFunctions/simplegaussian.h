#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction { //Gaussian thus extends WaveFunction
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double evaluate(std::vector<class Particle*> particles, int particle_id);
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles, int particle_id);
    //double computeDoubleDerivativeNumerically(std::vector<class Particle*> particles);
    std::vector<double> quantumForce(std::vector<class Particle*> particles);
    std::vector<double> quantumForce(std::vector<class Particle*> particles, int particle_id);
    void sample(std::vector<class Particle*> particles, double localEnergy);
    double totalRadius(std::vector<class Particle*> particles);
    double localEnergyTotalRadius(std::vector<class Particle*> particles, double localEnergy);
    void computeAverages(double steps);
    void gradientDescent();

protected:
  double m_av_total_radius = 0;
  double m_av_local_energy_total_radius = 0;
};
