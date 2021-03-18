#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToFile();
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    void setWriteout(bool samplertype){m_samplertype=samplertype;}
    void initiateFile                 ();
    void writeExpectationEnergyToFile (double cumul_energy, double local_energy);
    void writeLocalEnergyToFile      ();
    void setFileNameforEnergy(std:: string filename);

    virtual void gdsampler(std::vector<class Particle*> particles, double localEnergy){
      (void) particles;
      (void) localEnergy;
      return;
    }

protected:
    bool    m_samplertype=true;
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_kineticenergy = 0, m_cumulkinetic=0;
    double  m_potentialenergy = 0, m_cumulpotential=0;
    int     m_accepted=0;
    int     m_writeOutStep=100;
    class System* m_system = nullptr;
    std::string  m_energyfile;
    std::ofstream myfile;

};
