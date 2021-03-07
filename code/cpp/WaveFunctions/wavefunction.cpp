#include "wavefunction.h"

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setTolerance(double tol){
  m_tol = tol;
}
