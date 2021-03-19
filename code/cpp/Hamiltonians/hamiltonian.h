#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles)          = 0;//Compute the local energy
    virtual double computeLocalPotentialEnergy(std::vector<class Particle*> particles) = 0;//Compute the potential energy
    virtual double computeLocalKineticEnergy(std::vector<class Particle*> particles)   = 0;//Compute the kinetic energy

protected:
    class System* m_system = nullptr;
};
