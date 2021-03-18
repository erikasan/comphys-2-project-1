#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "sampler.h"
#include "GDsampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/numericgaussian.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include <string>
#include <omp.h>

using namespace std;

int main(int argc, char *argv[]) {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles  = 10;
    int numberOfSteps      = (int) 1e6;
    double omega           = 25;           // Oscillator frequency.
    double alpha           = 0.5;          // Variational parameter.
    double stepLength      = 0.1;          // Metropolis step length.
    double equilibration   = 0.1;          // Amount of the total steps used
    int num_threads=4;
    if (argc>1){
      try{
          numberOfDimensions = atoi(argv[1]);
          numberOfParticles  = atoi(argv[2]);
          numberOfSteps      = atoi(argv[3]);
          omega              = atoi(argv[4]);
          alpha              = atof(argv[5]);
          stepLength         = atof(argv[6]);
          equilibration      = atof(argv[7]);
          num_threads=4;
          seed=atoi(argv[8]);

      }
      catch (int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
      }
    }
    omp_set_num_threads(num_threads);
    int total_energy=0;
    int total_N=num_threads*numberOfSteps*(1-equilibration);
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nproc = omp_get_num_threads();
      System* system = new System(seed);

      system->m_energyfile="main"+to_string(id);
      system->setSampler               (new Sampler(system));
      system->setOmega(omega);
      system->setHamiltonian           (new HarmonicOscillator(system, omega));
      system->setWaveFunction          (new SimpleGaussian(system, alpha));
      system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));
      system->setEquilibrationFraction (equilibration);
      system->setStepLength            (stepLength);
      system->runMetropolisSteps       (numberOfSteps,false);
    }

    return 0;
}
