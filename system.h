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
            double timestep,
            unsigned int numberOfEquilibrationSteps);

    std::unique_ptr<class Sampler> runMetropolisSteps(
            double timestep,
            unsigned int numberOfMetropolisSteps,
            unsigned int MaxVariations);

    double computeLocalEnergy();
    const std::vector<double>& getWaveFunctionParameters();
    double wfDerivative();

private:
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;

    std::unique_ptr<class Hamiltonian> m_hamiltonian;
    std::unique_ptr<class WaveFunction> m_waveFunction;
    std::unique_ptr<class MonteCarlo> m_solver;
    std::vector<std::unique_ptr<class Particle>> m_particles;
};

