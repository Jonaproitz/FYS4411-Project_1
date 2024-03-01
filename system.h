#pragma once

#include <memory>
#include <vector>


class System {
public:
    System(
            std::unique_ptr<class Hamiltonian> hamiltonian,
            std::unique_ptr<class WaveFunction> waveFunction,
            std::unique_ptr<class MonteCarlo> solver,
            std::vector<std::unique_ptr<class Particle>> particles);

    void runEquilibrationSteps(
            double stepLength,
            unsigned int numberOfEquilibrationSteps,
            double TimeStep);

    std::unique_ptr<class Sampler> runMetropolisSteps(
            double stepLength,
            double TimeStep,
            unsigned int numberOfMetropolisSteps,
            unsigned int MaxVariations);

    double computeLocalEnergy();
    const std::vector<double>& getWaveFunctionParameters();
    void adjustAlpha(double adjust);

private:
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;

    std::unique_ptr<class Hamiltonian> m_hamiltonian;
    std::unique_ptr<class WaveFunction> m_waveFunction;
    std::unique_ptr<class MonteCarlo> m_solver;
    std::vector<std::unique_ptr<class Particle>> m_particles;
};

