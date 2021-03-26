#include "system.h"
#include <cassert>
#include "sampler.h"
#include "GDsampler.h"
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


bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test.
     */

     //Pick random Particle
     int particle_id = m_random->nextInt(m_numberOfParticles-1);

     //keep old position & wf
     double wf_old = m_waveFunction->evaluate(m_particles,particle_id);
     std::vector<double> position = m_particles[particle_id]->getPosition();

     //Move particle randomly
     for (int i = 0; i < m_numberOfDimensions; i++){
       m_particles[particle_id]->adjustPosition(2*(m_random->nextDouble() - 0.5)*m_stepLength,i);
     }
     m_waveFunction->updateDistances(m_particles,particle_id);
     double wf_new = m_waveFunction->evaluate(m_particles,particle_id);

     //Perform Metropolis test
     if ((wf_new*wf_new)/(wf_old*wf_old) > m_random->nextDouble() ){
       return true;
     }


     //Move particle back if metropolis test failed
     else{
       m_particles[particle_id]->setPosition(position);
       m_waveFunction->updateDistances(m_particles,particle_id);
     }
    return false;
}


void System::runMetropolisSteps(int numberOfMetropolisSteps, bool desire_output) {
    //Initiate information
    m_particles                 = m_initialState->getParticles();
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_waveFunction->initiateDistances(m_particles);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    //actually start the timer
    auto start = high_resolution_clock::now();

    bool acceptedStep;

    //warmup
    for (int i = 0; i < m_equilibrationSteps; i++) {
        acceptedStep = metropolisStep();
    }

    //run metropolis steps
    for (int i = 0; i < numberOfMetropolisSteps; i++) {
        acceptedStep = metropolisStep();
        m_sampler->sample(acceptedStep);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    m_duration = duration.count();

    m_sampler->computeAverages();
    if(m_sampler->getSamplePosition()){
      m_sampler->printPositions();
    }
    if (desire_output){
      m_sampler->printOutputToTerminal();
      m_sampler->printOutputToFile();
    }
}

void System::stopGradientDescent(){
  m_stopGradientDescent = true;
}

void System::gradientDescent(double tol, double learningRate, int maxIter){
  setSampler(new GDsampler(this));
  m_sampler -> setWriteout(false);
  m_tol = tol;
  m_waveFunction->setTolerance(tol);
  m_waveFunction->setLearningRate(learningRate);

  int fewMetropolisSteps = m_numberOfMetropolisSteps;

  int i = 0;
  while ((m_stopGradientDescent == false) & (i < maxIter)){
    runMetropolisSteps(fewMetropolisSteps, false);
    m_waveFunction->gradientDescent();
    i++;
  }

  runMetropolisSteps(m_numberOfMetropolisSteps, false); // All gas no brakes

  return;
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
    m_stepLengthRoot = sqrt(stepLength);
}

void System::setEquilibrationSteps(int equilibrationSteps) {
    assert(equilibrationSteps >= 0);
    m_equilibrationSteps = equilibrationSteps;
}

void System::setOmega(double omega){
  m_omega=omega;
}

void System::setMetropolisSteps(int numberOfMetropolisSteps){
  m_numberOfMetropolisSteps = numberOfMetropolisSteps;
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

void System::setSampler(Sampler* sampler){
  m_sampler = sampler;
}
