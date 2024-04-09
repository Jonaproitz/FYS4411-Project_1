#pragma once
#include <memory>
#include <vector>

class Hamiltonian {
public:
    // Use the default deconstructor generated by the compiler
    virtual ~Hamiltonian() = default;
    virtual double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
    ) = 0;
    virtual double computeOnebodyDensity(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
    ) = 0;
};

