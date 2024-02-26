#pragma once
#include <memory>
#include <vector>

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps);


    void sample(unsigned int acceptedStep, class System* system);
    void printOutputToTerminal(class System& system);
    void computeAverages();
    double getEnergy() { return m_energy.back(); }

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    std::vector<double> m_energy;
    std::vector<double> m_energy2;
    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergy2 = 0;
    double m_stepLength = 0;
};
