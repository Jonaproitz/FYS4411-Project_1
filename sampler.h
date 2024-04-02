#pragma once
#include <memory>
#include <vector>

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double timestep,
        unsigned int numberOfMetropolisSteps);


    void sample(unsigned int acceptedStep, class System* system);
    void printOutputToTerminal(class System& system);
    void printOutputToFile();
    void computeAverages();
    double getEnergy() { return m_energy; }
    double getVariance() {return m_variance; }
    double getCumulativeEnergy() {return m_cumulativeEnergy; }
    double getStepNumber() {return m_stepNumber; }

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    double m_energy = 0;
    double m_variance = 0;
    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergy2 = 0;
    double m_stepLength = 0;
};
