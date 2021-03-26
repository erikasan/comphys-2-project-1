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
#include <iomanip>
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
void Sampler::sample_position(){
  /*Sample the position and categorize it in increments of 0.05 for
  x between -5 and 5, radius between 0 and 10*/
  std::vector<Particle*> particles = m_system->getParticles();
  int numberOfParticles=m_system->getNumberOfParticles();
  double radius=0;
  double xpos=0;
  for (int i=0; i<numberOfParticles;i++){
    xpos=particles[i]->getPosition()[0];
    radius=sqrt(particles[i]->getRadiussquared()); //Oh the irony
    if(radius<10){
      rho_values[(int)lrint((radius*20))]++;
    }
    if((xpos<5) && (xpos>(-5))) {
      x_values[(int)lrint((xpos*20)+100)]++;
    }
  }
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

    m_stepNumber++;
    local_energyarray[m_stepNumber%m_writeOutStep]=localEnergy;
    if(m_stepNumber%m_writeOutStep==0 && m_samplertype==true){
      writeExpectationEnergyToFile();
    }
    if(m_sampleposition == true){
      sample_position();
    }
    m_accepted+=int(acceptedStep);
}
void Sampler::initiateFile                 (){
  myfile.open(m_energyfile,std::ios::out);
  myfile << "cumulative energy, local energy at each "<<m_writeOutStep << "step\n";
  myfile.close();
}
void Sampler::writeExpectationEnergyToFile (){
  //Write the whole array to file (every m_writeOutStep time)
  myfile.open(m_energyfile,std::ios::app);
  myfile << std::setprecision(10);
  for(int i=0;i<m_writeOutStep;i++){
    myfile <<local_energyarray[i] <<"\n";
  }
  myfile.close();
}
void Sampler::writeExpectationEnergyToFile (double local_energy){
  //Write the expectation local energy to file (one value at a time)
  myfile.open(m_energyfile,std::ios::app);
  myfile << std::setprecision(10) <<local_energy <<"\n";
  myfile.close();
}
void Sampler::setFileNameforEnergy(std::string filename){
  //Set filename for position and energy files for later analysis
  m_energyfile = m_system->getPath()+"energies_"+ filename+".csv";
  m_positionfile = m_system->getPath()+"positions_"+ filename+".csv";
}
void Sampler::printOutputToTerminal() {
    //Print the output to terminal
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
    cout << "standard error: " << sqrt(m_energysquared-m_energy*m_energy)/sqrt(m_stepNumber_count)<<endl;
    cout << "Duration: " << dur << " ms" << endl;
    cout << endl;
}
void Sampler::printOutputToFile(){
  //Prints output to file
  int np  = m_system->getNumberOfParticles();
  int nd  = m_system->getNumberOfDimensions();
  int ms  = m_system->getNumberOfMetropolisSteps();
  int p   = m_system->getWaveFunction()->getNumberOfParameters();
  (void)p;
  int dur = m_system->getDuration()/(1000);
  double ef = m_system->getEquilibrationSteps();
  std::vector<double> pa = m_system->getWaveFunction()->getParameters();
  std::ofstream newmyfile;
  std::ifstream infile(m_system->getPath()+"sympleharmonic.csv");
  if(infile.good()==false){
    newmyfile.open(m_system->getPath()+"sympleharmonic.csv",std::ofstream::app);
    newmyfile<<"Type,num_part,num_dim,steps,warmup,alpha,energy,kin_en,pot_en,acceptance_rate,time\n";
    newmyfile.close();
  }
  newmyfile.open(m_system->getPath()+"sympleharmonic.csv",std::ofstream::app);
  newmyfile << std::setprecision(10);
  newmyfile << "MonteCarlo,"<<np<<","<<nd<<","<<ms<<","<<ef;
  newmyfile <<","<<pa.at(0)<<","<<m_energy<<","<<m_kineticenergy<<","<<m_potentialenergy<<","<<double(m_accepted)/(ms)<<","<<dur<<endl;
  newmyfile.close();
  return;
}
void Sampler::printPositions() {
  //Unused in the final program, but possibly useful for later
  myfile.open(m_positionfile,std::ios::out);
  myfile << "xval,magnitude,radval,magnitude"<<endl;
  for(int i=0;i<rho_length;i++){
    myfile << -5+0.05*i<<","<<x_values[i] << ","<< 0.05*i<< ","<< rho_values[i]<<endl;
  }
  myfile.close();
}
void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities.*/
    double steps = m_stepNumber;
    m_energy = m_cumulativeEnergy/steps;
    m_kineticenergy = m_cumulkinetic/steps;
    m_potentialenergy = m_cumulpotential/steps;
    m_energysquared = m_cumulenergysquared/steps;
    m_system->getWaveFunction()->computeAverages(steps);
    m_stepNumber_count=m_stepNumber; //For presenting (integer)
    m_stepNumber = 0; //Set to zero for gradient descent (we start counting from 0)
    return;
}
