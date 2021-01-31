#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;

    /* This function, if not overridden in a subclass, evaluates the
    whole wave function (which makes no difference in terms of result
    (only timing) for Metropolis)*/
    virtual double evaluate(std::vector<class Particle*> particles, int particle_id) {(void)particle_id;return evaluate(particles);};
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};
