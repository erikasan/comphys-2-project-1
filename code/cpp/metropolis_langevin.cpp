#include "metropolis_langevin.h"
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



bool MetropolisLangevin::metropolisStep(){
  //Pick random Particle
  int particle_id = m_random->nextInt(m_numberOfParticles-1);
  //Evaluate old wave function, keep old position, evaluate quantum force
  double wf_old = m_waveFunction->evaluate(m_particles,particle_id);
  std::vector<double> position = m_particles[particle_id]->getPosition(); //old position
  std::vector<double> quantumForceOld = m_waveFunction->quantumForce(m_particles,particle_id);
  double adjustment = 0;
  //For each dimension, move particle
  for (int i = 0; i < m_numberOfDimensions; i++){
    adjustment = 0.5*quantumForceOld[i]*m_stepLength + m_random->nextGaussian(0,1)*m_stepLengthRoot;
    m_particles[particle_id]->adjustPosition(adjustment,i);
  }
  //Update distances to other particles (For SimpleGaussian, this does nothing)
  m_waveFunction->updateDistances(m_particles,particle_id);
  //Update position, evaluate new wave function and quantum force
  std::vector<double> newPosition = m_particles[particle_id]->getPosition();
  double wf_new = m_waveFunction->evaluate(m_particles,particle_id);
  std::vector<double> quantumForceNew = m_waveFunction->quantumForce(m_particles,particle_id);

  //Evaluate greens function
  double green = 0;
  for (int j = 0; j < m_numberOfDimensions; j++){
    green += 0.5*(quantumForceOld[j] + quantumForceNew[j])*
            (0.5*m_stepLength*0.5*(quantumForceOld[j] - quantumForceNew[j]) - newPosition[j] + position[j]);
  }

  //Metropolis-Hastings step
  if ((wf_new*wf_new)/(wf_old*wf_old)*exp(green) > m_random->nextDouble() ){
    return true;
  }

  else{
    m_particles[particle_id]->setPosition(position);
    m_waveFunction->updateDistances(m_particles,particle_id);
  }

  return false;
}

void MetropolisLangevin::runMetropolisSteps(int numberOfMetropolisSteps, bool desire_output){
  m_particles               = m_initialState->getParticles();
  m_numberOfMetropolisSteps = numberOfMetropolisSteps;

  //Initiate stuff
  m_waveFunction->initiateDistances(m_particles);
  m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
  //Start timer

  auto start = high_resolution_clock::now();
  bool acceptedStep;

  // Run m_equilibrationSteps equilibration steps
  for (int i = 0; i < m_equilibrationSteps; i++) {
      acceptedStep = metropolisStep();
  }
  //Run numberOfMetropolisSteps actual sampling steps
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
