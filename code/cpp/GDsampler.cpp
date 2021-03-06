#include "GDsampler.h"
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"

void GDsampler::gdsampler(std::vector<class Particle*> particles, double localEnergy){
  particles = m_system->getParticles();
  m_system->getWaveFunction()->sample(particles, localEnergy);
  return;
}
