#pragma once
#include "wavefunction.h"

class ComplexFunction : public WaveFunction { //ComplexFunction thus extends WaveFunction
public:
    ComplexFunction(class System* system, double alpha, double beta_param, double a_param);
    void initiateDistances(std::vector<class Particle*> particles);
    void updateDistances(std::vector<class Particle*> particles, int particle_id);
    double evaluate(std::vector<class Particle*> particles); //Evaluate the wave function part of all particles
    double evaluate(std::vector<class Particle*> particles, int particle_id); //Evaluate the wave function part of a fixed particle, requires seperability
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles, int particle_id);
    std::vector<double> quantumForce(std::vector<class Particle*> particles);
    std::vector<double> quantumForce(std::vector<class Particle*> particles, int particle_id);

    //
    void sample(std::vector<class Particle*> particles, double localEnergy);
    double totalRadius(std::vector<class Particle*> particles);
    double localEnergyTotalRadius(std::vector<class Particle*> particles, double localEnergy);
    void computeAverages(double steps);
    void gradientDescent();
    //

private:
    double beta=2.82843;
    double a=0.0043;
    //double beta=1;
    //double a=0.0000000000000001;

    double m_av_total_radius = 0;
    double m_av_local_energy_total_radius = 0;

    double ** particle_distances_absolute;

    double ** createNNMatrix(int n);
    void calculateParticleDistances(std::vector<class Particle*> particles); //Update for all particles
    void updateParticleDistances(std::vector<class Particle*> particles,int particle_id); //Update only for a given particle
    std::vector<double> distance(std::vector<class Particle*> particles, int particle_id1, int particle_id2);
    double distance(std::vector<double> part1, std::vector<double> part2);
    //std::vector<double> distance_vector(std::vector<double> part1, std::vector<double> part2);
    void printMatrix(double **A, int n);
    double u(double r);
    double uder(double r);
    double uderder(double r);

};
