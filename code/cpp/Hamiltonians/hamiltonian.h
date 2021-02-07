#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles)          = 0;
    virtual double computeLocalPotentialEnergy(std::vector<class Particle*> particles) = 0;
    virtual double computeLocalKineticEnergy(std::vector<class Particle*> particles)   = 0;

protected:
    class System* m_system = nullptr;
};
