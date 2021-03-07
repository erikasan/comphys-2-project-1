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

using namespace std;

int main(int argc, char *argv[]) {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles  = 1;
    int numberOfSteps      = (int) 1e6;
    double omega           = 25;           // Oscillator frequency.
    double alpha           = 1;          // Variational parameter.
    double stepLength      = 0.1;          // Metropolis step length.
    double equilibration   = 0.1;          // Amount of the total steps used

    if (argc>1){
      try{
          numberOfDimensions = atoi(argv[1]);
          numberOfParticles  = atoi(argv[2]);
          numberOfSteps      = atoi(argv[3]);
          omega              = atoi(argv[4]);
          alpha              = atof(argv[5]);
          stepLength         = atof(argv[6]);
          equilibration      = atof(argv[7]);
      }
      catch (int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
      }
    }

    System* system = new System(seed);
    system->setOmega(omega);
    system->setHamiltonian           (new HarmonicOscillator(system, omega));
    system->setWaveFunction          (new SimpleGaussian(system, alpha));
    system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setSampler               (new GDsampler(system));
    system->setEquilibrationFraction (equilibration);
    system->setStepLength            (stepLength);
    system->setMetropolisSteps       (numberOfSteps);

    double tol = 0.0001;
    int maxIter = 300;
    double learningRate = 0.01;
    system->gradientDescent(tol, learningRate, maxIter);

    cout << "Success!" << endl;
    return 0;
}
