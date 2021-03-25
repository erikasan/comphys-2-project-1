#include <iostream>
#include <string>
#include "../system.h"
#include "../particle.h"
#include "../metropolis_langevin.h"
#include "../sampler.h"
#include "../GDsampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/numericgaussian.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include "../WaveFunctions/complexfunction.h"
#include "../Hamiltonians/ellipticoscillator.h"
#include "../InitialStates/randomuniform2.h"
#include <omp.h>
#include <cmath>
using namespace std;

int main(int argc, char *argv[]) {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions = 3;
    int numberOfParticles  = 10;
    int numberOfSteps      = (int) 1e6;
    int num_threads=4;
    double alpha           = 0.5;          // Variational parameter.
    double beta             = 2.82843;
    double a                = 0.0043;
    double stepLength      = 0.1;          // Metropolis step length.
    int equilibration      = (int) 1e5;          // Amount of the total steps used
    string wF_type         = "HO"; // HO for Harmonic Oscillator, EO for Elliptic Oscillator, "NHO" for numerical harmonic oscillator
    string sampler_type    = "VMC"; //VMC for Brute Force, IMP for Importance sampling
    string filename_blocking = "no"; //no for "don't write", any other will then be written to file
    string path= "../../../output/";
    if (argc>=12){
          numberOfDimensions = atoi(argv[1]);
          numberOfParticles  = atoi(argv[2]);
          numberOfSteps      = atoi(argv[3]);
          alpha              = atof(argv[4]);
          stepLength         = atof(argv[5]);
          equilibration      = atoi(argv[6]);
          seed               = atoi(argv[7]);
          wF_type            = argv[8];
          sampler_type       = argv[9];
          filename_blocking  = argv[10];
          num_threads        = atoi(argv[11]);
          if(argc>=13){
            path               = argv[12];
          }
    }
    omp_set_num_threads(num_threads);
    double total_cumulativeEnergysquared=0;
    double total_cumulativeEnergy=0;
    int total_N=num_threads*numberOfSteps;
    double cumulativeEnergies[num_threads]={};
    double cumulativeEnergiesSquared[num_threads]={};
    if(sampler_type.compare("IMP") != 0 &&sampler_type.compare("VMC") != 0){
      cout << "Not a legal sampler type" << endl;
      return 1;
    }
    if(wF_type.compare("HO")!=0 && wF_type.compare("EO")!=0 && wF_type.compare("NHO")!=0){
        cout << "Not a legal Oscillator type" << endl;
    }
    #pragma omp parallel
    {
      System* system;
      int id = omp_get_thread_num();
      string filename_blocking_parallel=filename_blocking+to_string(id);
      int seed_parallel=seed+id;

      if (sampler_type.compare("VMC")==0){
        system = new System(seed_parallel);
      }
      else if (sampler_type.compare("IMP")==0){
        system = new MetropolisLangevin(seed_parallel);
      }

      system->m_energyfile=filename_blocking_parallel;
      system->setPath(path); //Change this is you want the output somewhere else
      system->setSampler               (new Sampler(system));
      system->setOmega(1);
      if(wF_type.compare("HO")==0){
        system->setHamiltonian           (new HarmonicOscillator(system, 1));
        system->setWaveFunction          (new SimpleGaussian(system, alpha));
        system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));
      }
      else if(wF_type.compare("NHO")==0){
        system->setHamiltonian           (new HarmonicOscillator(system, 1));
        system->setWaveFunction          (new NumericGaussian(system, alpha));
        system->setInitialState          (new RandomUniform(system, numberOfDimensions, numberOfParticles));
      }
      else if(wF_type.compare("EO")==0){
        system->setHamiltonian              (new EllipticOscillator(system, 1));
        system->setWaveFunction             (new ComplexFunction(system, alpha, beta, a));
        system->setInitialState             (new RandomUniformMinDist(system, numberOfDimensions, numberOfParticles,a));
      }

      system->setEquilibrationSteps (equilibration);
      system->setStepLength            (stepLength);
      if(filename_blocking.compare("no")==0){
        system->getSampler()->setWriteout(false);
      }
      system->runMetropolisSteps       (numberOfSteps,true);
      cumulativeEnergies[id]=system->getSampler()->getCumulEnergy();
      cumulativeEnergiesSquared[id]=system->getSampler()->getCumulEnergysquared();
    }
    for (int i=0; i<num_threads;i++){
      total_cumulativeEnergy+=cumulativeEnergies[i];
      total_cumulativeEnergysquared+=cumulativeEnergiesSquared[i];
    }
    double expectation_energy=total_cumulativeEnergy/total_N;
    double expectation_energysquared=total_cumulativeEnergysquared/total_N;
    cout << "Number of total steps:"<<total_N << endl;
    cout << "Energy: "<<expectation_energy<<endl;
    cout << "Energy squared: " <<expectation_energysquared <<endl;
    double std=sqrt(expectation_energysquared-expectation_energy*expectation_energy);
    cout << "standard Deviation: "<<std<<endl;
    cout << "standard error: " << std/sqrt(total_N)<<endl;
    return 0;
}
