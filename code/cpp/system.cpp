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

using namespace std::chrono;
System::System() {
    m_random = new Random(); // Initialize System with a random value
}

System::System(int seed) {
    m_random = new Random(seed); // Initialize System with a seed
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
     double wf_new=m_waveFunction->evaluate(m_particles,particle_id);
     //Change position randomly
     if ((wf_new*wf_new)/(wf_old*wf_old)> m_random->nextDouble() ){
       return true;
     }//Perform Metropolis test
     else{
       m_particles[particle_id]->setPosition(position);
     }
    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
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
    m_sampler->printOutputToTerminal();
    m_sampler->printOutputToFile();
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
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
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
