#pragma once
#include <vector>
#include "Math/random.h"
#include <string>
class System {
public:
    System();
    System(int seed);

    virtual bool metropolisStep       ();
    virtual void runMetropolisSteps   (int numberOfMetropolisSteps, bool desire_output);

    void gradientDescent              (double tol, double learningRate, int maxIter);

    void setNumberOfParticles         (int numberOfParticles);
    void setNumberOfDimensions        (int numberOfDimensions);

    void setStepLength                (double stepLength);
    void setEquilibrationSteps        (int equilibrationSteps);

    void setHamiltonian               (class Hamiltonian* hamiltonian);
    void setWaveFunction              (class WaveFunction* waveFunction);
    void setInitialState              (class InitialState* initialState);
    void setSampler                   (class Sampler* sampler);

    void setOmega                     (double omega);
    void setMetropolisSteps           (int numberOfMetropolisSteps);

    void stopGradientDescent          ();

    class WaveFunction* getWaveFunction(){
      return m_waveFunction;
    }

    class Hamiltonian* getHamiltonian(){
      return m_hamiltonian;
    }

    class Sampler* getSampler(){
      return m_sampler;
    }

    std::vector<class Particle*> getParticles(){
      return m_particles;
    }

    class Random* getRandomEngine(){
      return m_random;
    }

    int getNumberOfParticles(){
      return m_numberOfParticles;
    }

    int getNumberOfDimensions(){
      return m_numberOfDimensions;
    }

    int getNumberOfMetropolisSteps(){
      return m_numberOfMetropolisSteps;
    }

    int getDuration(){
      return m_duration;
    }

    int getEquilibrationSteps(){
      return m_equilibrationSteps;
    }

    double getOmega(){
      return m_omega;
    }
    std::string getPath(){
      return m_path;
    }
    void setPath(std::string path){
      m_path=path;
    }
    std::string  m_energyfile="default";
    std::string  m_path="../../../output/";
protected:
    int                          m_numberOfParticles       = 0;
    int                          m_numberOfDimensions      = 0;
    int                          m_numberOfMetropolisSteps = 0;
    int                          m_equilibrationSteps      = 0.0;
    double                       m_stepLength              = 0.1;
    double                       m_stepLengthRoot          = 0.01;
    int                          m_duration                = 0;
    double                       m_omega                   = 0.0;
    double                       m_tol                     = 0.0;
    bool                         m_stopGradientDescent     = false;

    class WaveFunction*          m_waveFunction = nullptr;
    class Hamiltonian*           m_hamiltonian  = nullptr;
    class InitialState*          m_initialState = nullptr;
    class Sampler*               m_sampler      = nullptr;
    std::vector<class Particle*> m_particles    = std::vector<class Particle*>();
    class Random*                m_random       = nullptr;
};
