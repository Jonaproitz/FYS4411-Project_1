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
    void computeAverages();
    double getEnergy() { return m_energy.back(); }
    double getVariance() {return m_variance.back(); }
    void printOutputToFile();
    void storeAlphaValues(double a) {m_alphaValues.push_back(a);}

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    std::vector<double> m_energy;
    std::vector<double> m_variance;
    std::vector<double> m_alphaValues;
    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergy2 = 0;
    double m_stepLength = 0;
};
