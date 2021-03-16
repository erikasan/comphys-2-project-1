#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;

    /* This function, if not overridden in a subclass, evaluates the
    whole wave function (which makes no difference in terms of result
    (only timing) for Metropolis)*/
    virtual double evaluate(std::vector<class Particle*> particles, int particle_id) {(void)particle_id;return evaluate(particles);};
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> quantumForce(std::vector<class Particle*> particles, int particle_id) {(void)particles;(void)particle_id;return quantumForce(particles);};
    virtual std::vector<double> quantumForce(std::vector<class Particle*> particles) = 0;
    virtual void updateDistances(std::vector<class Particle*> particles,int particle_id){(void)particle_id;(void)particles;return;}
    virtual void initiateDistances(std::vector<class Particle*> particles){(void)particles;return;}
    virtual void sample(std::vector<class Particle*> particles, double localEnergy){(void) localEnergy;(void) particles; return;};
    virtual void computeAverages(double steps){(void) steps; return;}
    virtual void gradientDescent(){return;}
    void setTolerance(double tol);
    void setLearningRate(double learningRate);


protected:
    int m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    double m_tol = 0.0;
    double m_learningRate = 0.0;
};
