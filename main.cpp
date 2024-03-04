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
    unsigned int numberOfDimensions = 1;
    unsigned int numberOfParticles = 1;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e5;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.2; // Variational parameter.
    double stepLength = 1; // Metropolis step length.

    // Set number of variations in alpha
    unsigned int MaxVariations = 70;

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(stepLength, numberOfDimensions, numberOfParticles, *rng);
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
            stepLength,
            numberOfEquilibrationSteps);

    // Run Metropolis algoritm
    auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps,
            MaxVariations);

    // Output information from the simulation to terminal
    sampler->printOutputToTerminal(*system);
    sampler->printOutputToFile();

    return 0;
}
