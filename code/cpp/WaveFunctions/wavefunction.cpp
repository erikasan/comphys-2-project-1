#include "wavefunction.h"

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setTolerance(double tol){
  m_tol = tol;
}

void WaveFunction::setLearningRate(double learningRate){
  m_learningRate = learningRate;
}
