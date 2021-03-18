#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../metropolis_langevin.h"
#include "sampler.h"
#include "GDsampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/complexfunction.h"
#include "../Hamiltonians/ellipticoscillator.h"
#include "../InitialStates/randomuniform2.h"

#include "../Hamiltonians/hamiltonian.h"
#include "../InitialStates/initialstate.h"
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
    double alpha           = 0.5;          // Variational parameter.
    double beta            = 2.82843;
    double a               = 0.0043;
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

    System* system = new MetropolisLangevin(seed);
    system->setOmega(omega);
    // system->setHamiltonian           (new HarmonicOscillator(system, omega));
    // system->setWaveFunction          (new SimpleGaussian(system, alpha));
    // system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    system->setHamiltonian              (new EllipticOscillator(system, omega));
    system->setWaveFunction             (new ComplexFunction(system, alpha, beta, a));
    system->setInitialState             (new RandomUniformMinDist(system, numberOfDimensions, numberOfParticles,a));


    system->setEquilibrationFraction (equilibration);
    system->setStepLength            (stepLength);
    system->setMetropolisSteps       (numberOfSteps);

    double tol = 0.0001;
    int maxIter = 100;
    double learningRate = 0.01;
    system->gradientDescent(tol, learningRate, maxIter);

    return 0;
}
