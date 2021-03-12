/*
compile as:
c++ -c unit_testing.cpp
c++ -o test.exe unit_testing.o VMC.o System.o functions.o tests_main.o
./test.exe
*/

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <iostream>
#include "../system.h"
#include "../metropolis_langevin.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/numericgaussian.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include "../sampler.h"
#include "../catch.hpp"
using namespace std;
TEST_CASE("Test that the analytical value matches the calculated value when the right omega is chosen"){
  int seed = 2020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 1e6;
  double omega            = 10;          // Oscillator frequency.
  double alpha            = 0.5;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  double equilibration    = 0.1;          // Amount of the total steps used
  string sample_type="not_numerically";
  System* system = new System(seed);
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new SimpleGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationFraction    (equilibration);
  system->setStepLength               (stepLength);
  system->runMetropolisSteps          (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*0.5;
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-3);
}
TEST_CASE("Evaluate wether numerical also works"){
  int seed = 2020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 1e6;
  double omega            = 10;          // Oscillator frequency.
  double alpha            = 0.5;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  double equilibration    = 0.1;          // Amount of the total steps used
  System* system = new System(seed);
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new NumericGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationFraction    (equilibration);
  system->setStepLength               (stepLength);
  system->runMetropolisSteps          (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*0.5;
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-3);
}
TEST_CASE("Evaluate wether importance sampling works"){
  int seed = 020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 1e5;
  double omega            = 20;          // Oscillator frequency.
  double alpha            = 0.5;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  double equilibration    = 0.1;          // Amount of the total steps used
  System* system = new MetropolisLangevin(seed);
  system->setOmega(omega);
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new SimpleGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationFraction    (equilibration);
  system->setStepLength               (stepLength);
  system->runMetropolisSteps  (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*0.5;
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-3);
}
