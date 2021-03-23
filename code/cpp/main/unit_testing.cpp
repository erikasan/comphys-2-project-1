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
#include "../WaveFunctions/complexfunction.h"
#include "../Hamiltonians/ellipticoscillator.h"
#include "../InitialStates/randomuniform2.h"
#include "../Math/random.h"
#include "../sampler.h"
#include "../GDsampler.h"
#include "../catch.hpp"
using namespace std;
TEST_CASE("Evaluate that Gradient Descent leads us to the correct minimum"){
  int seed = 2020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 5e4;
  double alpha            = 0.55;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  int equilibration      = (int) 1e4;          // Amount of the total steps used
  System* system = new MetropolisLangevin(seed);
  system->setSampler               (new Sampler(system));
  system->setOmega(25);
  system->setHamiltonian           (new HarmonicOscillator(system, 25));
  system->setWaveFunction          (new SimpleGaussian(system, alpha));
  system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationSteps (equilibration);
  system->setStepLength            (stepLength);
  system->setMetropolisSteps       (numberOfSteps);
  double tol = 0.00000001;
  int maxIter = 50;
  double learningRate = 0.02;
  system->gradientDescent(tol, learningRate, maxIter);
  double alpha_calculated=system->getWaveFunction()->getParameters()[0];
  double alpha_expected=0.5;
  REQUIRE(fabs(alpha_expected-alpha_calculated)<3e-3);
}
TEST_CASE("Test that the analytical value matches the calculated value when the right omega is chosen"){
  int seed = 2020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 1e6;
  double omega            = 10;          // Oscillator frequency.
  double alpha            = 1.0;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  int equilibration      = (int) 1e5;          // Amount of the total steps used
  System* system = new System(seed);
  system->setSampler               (new Sampler(system));
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new SimpleGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationSteps   (equilibration);
  system->setStepLength               (stepLength);
  system->getSampler()->setWriteout(false);
  system->runMetropolisSteps          (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*(alpha+(1-4*alpha*alpha)/(8*alpha));
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-1);
}

TEST_CASE("Evaluate wether numerical also works"){
  int seed = 2020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 10;
  int numberOfSteps       = (int) 1e6;
  double omega            = 10;          // Oscillator frequency.
  double alpha            = 1.0;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  int equilibration      = (int) 1e5;         // Amount of the total steps used
  System* system = new System(seed);
  system->setSampler               (new Sampler(system));
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new NumericGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationSteps    (equilibration);
  system->setStepLength               (stepLength);
  system->getSampler()->setWriteout(false);
  system->runMetropolisSteps          (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*(alpha+(1-4*alpha*alpha)/(8*alpha));
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-1);
}
TEST_CASE("Evaluate wether importance sampling works"){
  int seed = 020;

  int numberOfDimensions  = 3;
  int numberOfParticles   = 1;
  int numberOfSteps       = (int) 1e6;
  double omega            = 20;          // Oscillator frequency.
  double alpha            = 1.0;          // Variational parameter.
  double stepLength       = 0.1;          // Metropolis step length.
  int equilibration      = (int) 1e5;          // Amount of the total steps used
  System* system = new MetropolisLangevin(seed);
  system->setSampler                  (new Sampler(system));
  system->setOmega(omega);
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new SimpleGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationSteps    (equilibration);
  system->setStepLength               (stepLength);
  system->getSampler()->setWriteout(false);
  system->runMetropolisSteps  (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfDimensions*numberOfParticles*(alpha+(1-4*alpha*alpha)/(8*alpha));
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-1);
}
TEST_CASE("Test wether the difficult wave function, in the limit a=0, becomes a harmonic oscillator"){
  int seed = 2020;

  int numberOfDimensions = 3;
  int numberOfParticles  = 10;
  int numberOfSteps      = (int) 2e6;
  double alpha           = 1.0;          // Variational parameter.
  double beta             = 2.82843;
  double a                = 0; //teeny tiny number
  double stepLength      = 0.1;          // Metropolis step length.
  double multiplicator  = (1+1+beta);
  int equilibration      = (int) 1e5;          // Amount of the total steps used
  System* system = new MetropolisLangevin(seed);
  system->setSampler                  (new Sampler(system));
  system->setHamiltonian              (new EllipticOscillator(system, 1));
  system->setWaveFunction             (new ComplexFunction(system, alpha, beta, a));
  system->setInitialState             (new RandomUniformMinDist(system, numberOfDimensions, numberOfParticles,a));
  system->setEquilibrationSteps   (equilibration);
  system->setStepLength               (stepLength);
  system->getSampler()->setWriteout(false);
  system->runMetropolisSteps          (numberOfSteps,false);
  double potential_energy_calculated=system->getSampler()->getEnergy();
  double potential_energy_expected=numberOfParticles*multiplicator*(alpha+(1-4*alpha*alpha)/(8*alpha));
  REQUIRE(fabs(potential_energy_expected-potential_energy_calculated)<1e-1);
}
