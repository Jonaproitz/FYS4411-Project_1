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

void System::energyDerivative(
        double timestep,
        unsigned int numberOfMetropolisSteps)
{

    double energy = 0, deltaPsi = 0, derivativePsi = 0, energyDer = 0;
    double eta = 0.05, etol = 1e-4;
    unsigned int MaxVariations = 1000;
        for (unsigned int iter=0; iter<MaxVariations; iter++) {
            for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
                for (unsigned int j = 0; j < m_numberOfParticles; j++) {
                    for (unsigned int d = 0; d < m_numberOfDimensions; d++) {
                        // Call solver method to do a single Monte-Carlo step.
                        m_solver->step(timestep, *m_waveFunction, m_particles, j, d);
                    }
                }
                
                auto localEnergy = computeLocalEnergy();
                energy += localEnergy;

                auto WFder = wfDerivative();
                deltaPsi += WFder;
                derivativePsi += WFder*localEnergy;

            }
        energy /= numberOfMetropolisSteps;
        derivativePsi /= numberOfMetropolisSteps;
        deltaPsi /= numberOfMetropolisSteps;

        energyDer = 2*(derivativePsi-deltaPsi*energy);
        m_waveFunction->adjustAlpha(-eta*energyDer);
        if (abs(eta*energyDer) < etol) {
            std::cout << "iter = " << iter + 1 << std::endl;
            std::cout << "energyDer = " << energyDer << std::endl;
            break;}
    }
    std::cout << "minimal alpha: " << m_waveFunction->getParameters().at(0) << std::endl;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
        double timestep,
        unsigned int numberOfMetropolisSteps)
{
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions,
            timestep,
            numberOfMetropolisSteps);

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

double System::wfDerivative()
{
    // Helper function
    return m_waveFunction->wfDerivative(m_particles);
}
