#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2024;

    // Initial setup for simulation
    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 10;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e4;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e5;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.7; // Variational parameter.
    double timestep = 0.05; // Metropolis step length.

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(timestep, numberOfDimensions, numberOfParticles, *rng);
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            std::make_unique<HarmonicOscillator>(omega),
            // Construct unique_ptr to wave function
            std::make_unique<SimpleGaussian>(alpha),
            // Construct unique_ptr to solver, and move rng
            std::make_unique<Metropolis>(std::move(rng)),
            // Move the vector of particles to system
            std::move(particles));

    // Run steps to equilibrate particles
    system->runEquilibrationSteps(
            timestep,
            numberOfEquilibrationSteps);

    system->optimizeParameters(
            timestep,
            numberOfMetropolisSteps);

    numberOfMetropolisSteps = 1<<19;

    // Run Metropolis algoritm
    auto sampler = system->runMetropolisSteps(
            timestep,
            numberOfMetropolisSteps);

    // Output information from the simulation to terminal
    sampler->printOutputToTerminal(*system);


    return 0;
}
