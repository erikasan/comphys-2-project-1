#pragma once
#include <vector>

class Particle {
public:
    Particle();
    Particle(const Particle &particle);
    void setPosition(const std::vector<double> &position); //Function to set the particle's position
    void adjustPosition(double change, int dimension); //Function to update the position at dim dimension by a given change
    void setNumberOfDimensions(int numberOfDimensions); //Defines the number of dimensions the particle exists in
    std::vector<double> getPosition() { return m_position; } //Gets the particles Position
    double getRadiussquared(); //Get the Radius (distance from origin)

private:
    int     m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
};
