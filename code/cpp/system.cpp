#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>

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
     double wf_old=m_waveFunction->evaluate(m_particles);
     //Pick random Particle
     int particle_id = m_random -> nextInt(m_numberOfParticles-1);
     //keep old positions
     //std::cout << "test nr 1\n";
     std::vector<double> position= m_particles[particle_id]->getPosition();
     for (int i=0; i<m_numberOfDimensions;i++){
       m_particles[particle_id]->adjustPosition(2*(m_random->nextDouble()-0.5)*m_stepLength,i);
     }
     //std::cout << m_particles[particle_id]->getPosition()[0] << "\n";
     double wf_new=m_waveFunction->evaluate(m_particles);
     //Change position randomly
     if ((wf_new*wf_new)/(wf_old*wf_old)> m_random->nextDouble() ){
       //std::cout<< "WF_old: " << wf_old << "WF_new: " << wf_new <<"\n";
       return true;
     }//Perform Metropolis test
     else{
       //std::cout << "not accepted\n";
       //std::cout << "new X:" << m_particles[particle_id]->getPosition()[0]<<"\n";
       m_particles[particle_id]->setPosition(position);
       //std::cout << "old X:" << m_particles[particle_id]->getPosition()[0]<<"\n";
     }
    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

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
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
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
