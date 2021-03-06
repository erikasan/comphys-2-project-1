#pragma once
#include <vector>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToFile();
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }

    virtual void gdsampler(std::vector<class Particle*> particles, double localEnergy){
      (void) particles;
      (void) localEnergy;
      return;
    }

protected:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_kineticenergy = 0, m_cumulkinetic=0;
    double  m_potentialenergy = 0, m_cumulpotential=0;
    class System* m_system = nullptr;
};
