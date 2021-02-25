#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <chrono>
#include <string>
#include <cmath>
using namespace std::chrono;
using namespace std;
System::System() {
    m_random = new Random(); // Initialize System with a random value
}

System::System(int seed) {
    m_random = new Random(seed); // Initialize System with a seed
}
bool System::metropolis_LangevinStep(){
  //Pick random Particle
  int particle_id = m_random -> nextInt(m_numberOfParticles-1);
  double wf_old=m_waveFunction->evaluate(m_particles,particle_id);
  std::vector<double> position= m_particles[particle_id]->getPosition(); //old position
  std::vector<double> quantumForceOld= m_waveFunction->quantumForce(m_particles,particle_id);
  double adjustment=0;
  for (int i=0;i<m_numberOfDimensions;i++){
    adjustment=0.5*quantumForceOld[i]*m_stepLength+m_random->nextGaussian(0,1)*m_stepLengthRoot;
    m_particles[particle_id]->adjustPosition(adjustment,i);
  }
  m_waveFunction->updateDistances(m_particles,particle_id);
  std::vector<double> newPosition=m_particles[particle_id]->getPosition();
  double wf_new=m_waveFunction->evaluate(m_particles,particle_id);
  std::vector<double> quantumForceNew= m_waveFunction->quantumForce(m_particles,particle_id);
  double green=0;
  for (int j=0;j<m_numberOfDimensions;j++){
    green+=0.5*(quantumForceOld[j]+quantumForceNew[j])*
            (0.5*m_stepLength*0.5*(quantumForceOld[j]-quantumForceNew[j])-newPosition[j]+position[j]);
  }
  if ((wf_new*wf_new)/(wf_old*wf_old)*exp(green)> m_random->nextDouble() ){

    return true;
  }
  else{
      m_particles[particle_id]->setPosition(position);
      m_waveFunction->updateDistances(m_particles,particle_id);
  }
  return false;
}
bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

     //Pick random Particle
     int particle_id = m_random -> nextInt(m_numberOfParticles-1);
     double wf_old=m_waveFunction->evaluate(m_particles,particle_id);
     //keep old position
     std::vector<double> position= m_particles[particle_id]->getPosition();
     for (int i=0; i<m_numberOfDimensions;i++){
       m_particles[particle_id]->adjustPosition(2*(m_random->nextDouble()-0.5)*m_stepLength,i);
     }
     m_waveFunction->updateDistances(m_particles,particle_id);
     double wf_new=m_waveFunction->evaluate(m_particles,particle_id);
     //Change position randomly
     if ((wf_new*wf_new)/(wf_old*wf_old)> m_random->nextDouble() ){
       return true;
     }//Perform Metropolis test
     else{
       m_particles[particle_id]->setPosition(position);
       m_waveFunction->updateDistances(m_particles,particle_id);
     }
    return false;
}
void System::runMetropolisLangevinSteps(int numberOfMetropolisSteps, bool desire_output){
  m_particles                 = m_initialState->getParticles();
  m_sampler                   = new Sampler(this);
  m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
  m_waveFunction-> initiateDistances(m_particles);
  m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
  auto start = high_resolution_clock::now();
  for (int i=0; i < numberOfMetropolisSteps; i++) {
      bool acceptedStep = metropolis_LangevinStep();
      if (i > (int)(m_equilibrationFraction*numberOfMetropolisSteps)){
          m_sampler->sample(acceptedStep);
      }
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  m_duration=duration.count();
  m_sampler->computeAverages();
  if (desire_output){
    m_sampler->printOutputToTerminal();
    m_sampler->printOutputToFile();
  }
}
void System::runMetropolisSteps(int numberOfMetropolisSteps, bool desire_output) {
    m_particles                 = m_initialState->getParticles();

    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_waveFunction-> initiateDistances(m_particles);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    auto start = high_resolution_clock::now();
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        if (i > (int)(m_equilibrationFraction*numberOfMetropolisSteps)){
          m_sampler->sample(acceptedStep);
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    m_duration=duration.count();
    m_sampler->computeAverages();
    if (desire_output){
      m_sampler->printOutputToTerminal();
      m_sampler->printOutputToFile();
    }
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
    m_stepLengthRoot=sqrt(stepLength);
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}
void System::setOmega(double omega){
  m_omega=omega;
}
void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}
