#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void sample_numerically(bool acceptedStep);
    void printOutputToFile();
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_kineticenergy = 0, m_cumulkinetic=0;
    double  m_potentialenergy = 0, m_cumulpotential=0;
    class System* m_system = nullptr;
};
