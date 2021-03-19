#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <chrono>
#include <string>
using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
    setFileNameforEnergy(system->m_energyfile);
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {

    std::vector<Particle*> particles = m_system->getParticles();

    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulkinetic=0;
        m_cumulpotential=0;


        if (m_samplertype==true){
          initiateFile();
          }
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */

    double localPotentialEnergy = m_system->getHamiltonian()->
                         computeLocalPotentialEnergy(particles);
    double localKineticEnergy   = m_system->getHamiltonian()->
                         computeLocalKineticEnergy(particles);

    double localEnergy = localPotentialEnergy + localKineticEnergy;

    gdsampler(particles, localEnergy);
    m_cumulenergysquared   += localEnergy*localEnergy;
    m_cumulativeEnergy += localEnergy;
    m_cumulkinetic     += localKineticEnergy;
    m_cumulpotential   += localPotentialEnergy;
    if(m_stepNumber%m_writeOutStep==0 && m_samplertype==true){
      writeExpectationEnergyToFile(m_cumulativeEnergy/(m_stepNumber+1),localEnergy);
    }
    m_stepNumber++;
    m_accepted+=int(acceptedStep);
}
void Sampler::initiateFile                 (){
  myfile.open(m_energyfile,std::ios::out);
  myfile << "cumulative energy, local energy at each "<<m_writeOutStep << "step\n";
  myfile.close();
  myfile.open(m_energyfile,std::ios::app);
}
void Sampler::writeExpectationEnergyToFile (double cumul_energy, double local_energy){

  myfile <<local_energy <<"\n";
}
void Sampler::setFileNameforEnergy(std::string filename){
  m_energyfile = "../../../output/energies_"+ filename+".csv";
}
void Sampler::printOutputToTerminal() {
    int    np = m_system->getNumberOfParticles();
    int    nd = m_system->getNumberOfDimensions();
    int    ms = m_system->getNumberOfMetropolisSteps();
    int    p  = m_system->getWaveFunction()->getNumberOfParameters();
    int    dur= m_system->getDuration()/(1000);
    double ef = m_system->getEquilibrationSteps();

    std::vector<double> pa = m_system->getWaveFunction()->getParameters();
    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << "Acceptance rate:" << double(m_accepted)/(ms)<<  endl;
    cout << " Energy : " << m_energy << endl;
    cout << "Kinetic Energy: " << m_kineticenergy <<endl;
    cout << "Potential Energy: " << m_potentialenergy <<endl;
    cout << "standard error: " << sqrt(m_energysquared-m_energy*m_energy)/sqrt(m_stepNumber)<<endl;
    cout << "Duration: " << dur << " ms" << endl;
    cout << endl;
}
void Sampler::printOutputToFile(){
  int np  = m_system->getNumberOfParticles();
  int nd  = m_system->getNumberOfDimensions();
  int ms  = m_system->getNumberOfMetropolisSteps();
  int p   = m_system->getWaveFunction()->getNumberOfParameters();
  int dur = m_system->getDuration()/(1000);
  double ef = m_system->getEquilibrationSteps();
  std::vector<double> pa = m_system->getWaveFunction()->getParameters();
  std::ofstream myfile;
  myfile.open("../../../output/sympleharmonic.csv",std::ofstream::app);
  myfile << "sympleharmonic,"<<np<<","<<nd<<","<<ms<<","<<ef<<",";
  myfile <<pa.at(0)<<","<<","<<m_energy<<","<<m_kineticenergy<<","<<m_potentialenergy<<","<<double(m_accepted)/(ms)<<","<<dur<<endl;
  myfile.close();
  std::cout << "written to file...?"<<endl;
  return;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities.*/
    double steps = m_stepNumber;
    m_energy = m_cumulativeEnergy/steps;
    m_kineticenergy = m_cumulkinetic/steps;
    m_potentialenergy = m_cumulpotential/steps;
    m_energysquared = m_cumulenergysquared/steps;
    m_system->getWaveFunction()->computeAverages(steps);
    cout << m_energy << endl;
    myfile.close();
    return;
}
