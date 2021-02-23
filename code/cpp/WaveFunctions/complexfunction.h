#pragma once
#include "wavefunction.h"

class ComplexFunction : public WaveFunction { //Gaussian thus extends WaveFunction
public:
    ComplexFunction(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double evaluate(std::vector<class Particle*> particles, int particle_id);
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles, int particle_id);
    //double computeDoubleDerivativeNumerically(std::vector<class Particle*> particles);
    std::vector<double> quantumForce(std::vector<class Particle*> particles);
    std::vector<double> quantumForce(std::vector<class Particle*> particles, int particle_id);

private:
    double beta=2.82843;
    double a=0.0043;
    double ** particle_distances_absolute;

    double ** createNNMatrix(int n);
    void calculateParticleDistances(std::vector<class Particle*> particles); //Update for all particles
    void updateParticleDistances(std::vector<class Particle*> particles,int particle_id); //Update only for a given particle
    std::vector<double> distance(std::vector<class Particle*> particles, int particle_id1, int particle_id2);
    double vectorProduct(std::vector<double> part1, std::vector<double> part2);
    double distance(std::vector<double> part1, std::vector<double> part2);
    std::vector<double> distance_vector(std::vector<double> part1, std::vector<double> part2);
};
