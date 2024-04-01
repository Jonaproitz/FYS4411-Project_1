#include <iostream>
#include <memory>
#include <cassert>

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"


System::System(
        std::unique_ptr<class Hamiltonian> hamiltonian,
        std::unique_ptr<class WaveFunction> waveFunction,
        std::unique_ptr<class MonteCarlo> solver,
        std::vector<std::unique_ptr<class Particle>> particles)
{
    m_numberOfParticles = particles.size();
    m_numberOfDimensions = particles[0]->getNumberOfDimensions();
    m_hamiltonian = std::move(hamiltonian);
    m_waveFunction = std::move(waveFunction);
    m_solver = std::move(solver);
    m_particles = std::move(particles);
}


void System::runEquilibrationSteps(
        double timeStep,
        unsigned int numberOfEquilibrationSteps)
{
    for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++) {
        for (unsigned int j = 0; j < m_numberOfParticles; j++) {
            for (unsigned int d = 0; d < m_numberOfDimensions; d++) {
                m_solver->step(timeStep, *m_waveFunction, m_particles, j, d);
            }
        }
    }

    return;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
        double timestep,
        unsigned int numberOfMetropolisSteps,
        unsigned int MaxVariations)
{
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions,
            timestep,
            numberOfMetropolisSteps);
    // Set adjustment to alpha
    double adjust = 0.02;
    for (unsigned int m = 0; m <= MaxVariations; m++) {
        // Store alpha value
        sampler->storeAlphaValues(getWaveFunctionParameters().at(0));
        for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
            unsigned int numberOfAcceptedSteps = 0;
            for (unsigned int j = 0; j < m_numberOfParticles; j++) {
                for (unsigned int d = 0; d < m_numberOfDimensions; d++) {
                    // Call solver method to do a single Monte-Carlo step.
                    bool acceptedStep = m_solver->step(timestep, *m_waveFunction, m_particles, j, d);
                    numberOfAcceptedSteps += acceptedStep;
                }
            }
        
        // Sample the energy
        sampler->sample(numberOfAcceptedSteps, this);
    }

    sampler->computeAverages();
    adjustAlpha(adjust);
    }
    return sampler;
}

double System::computeLocalEnergy()
{
    // Helper function
    return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}

const std::vector<double>& System::getWaveFunctionParameters()
{
    // Helper function
    return m_waveFunction->getParameters();
}

void System::adjustAlpha(double adjust)
{
    // Helper function
    m_waveFunction->adjustAlpha(adjust);
}