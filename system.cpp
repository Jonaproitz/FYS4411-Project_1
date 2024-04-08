#include <iostream>
#include <memory>
#include <cassert>
#include <fstream>

#include "/opt/homebrew/opt/libomp/include/omp.h"

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"

#include <iomanip>


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
            // Call the solver method to do a single Monte Carlo step for each particle
            m_solver->step(timeStep, *m_waveFunction, m_particles, j);
        }
    }

    return;
}

void System::optimizeParameters(
        double timestep,
        unsigned int numberOfMetropolisSteps)
{
    // Set endpoint and weighting of optimization
    double etol = 1e-2, eta = 1e-2 / m_numberOfParticles;
    unsigned int MaxVariations = 20;
    for (unsigned int iter=0; iter<MaxVariations; iter++) {
        // Define properties used in optimization
        double energy = 0, deltaPsi = 0, derivativePsi = 0, energyDer = 0;
        for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
            for (unsigned int j = 0; j < m_numberOfParticles; j++) {
                // Call solver method to do a single Monte-Carlo step.
                m_solver->step(timestep, *m_waveFunction, m_particles, j);
            }
            // Sample necessary properties to find the energy derivative
            auto localEnergy = computeLocalEnergy();
            energy += localEnergy;

            auto WFder = wfDerivative();
            deltaPsi += WFder;
            derivativePsi += WFder*localEnergy;
        }
        // Average results of the run
        energy /= numberOfMetropolisSteps;
        derivativePsi /= numberOfMetropolisSteps;
        deltaPsi /= numberOfMetropolisSteps;
        
        // Find Energy derivative
        energyDer = 2*(derivativePsi-deltaPsi*energy);
        const char separator    = ' ';
        const int nameWidth     = 20;
        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << m_waveFunction->getParameters().at(0);
        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << energy/m_particles.size();
        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << energyDer << std::endl;
        
        if (abs(energyDer) < etol) {
            std::cout << "iter = " << iter + 1 << std::endl;
            std::cout << "energyDer = " << energyDer << std::endl;
            break;}

        // Adjust alpha value
        double adjust = -energyDer*eta;
        m_waveFunction->adjustAlpha(adjust);
    }
    std::cout << "minimal alpha: " << m_waveFunction->getParameters().at(0) << std::endl;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
        double timestep,
        unsigned int numberOfMetropolisSteps)
{
    unsigned int numberOfThreads = 4;
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions,
            timestep,
            numberOfMetropolisSteps*numberOfThreads);
    
    std::ofstream myfile;
    myfile.open("Results/Energies_G_" + std::to_string(m_numberOfParticles) + ".dat");
    
    // Start parallelization
    #pragma omp parallel num_threads(numberOfThreads)
    {
    unsigned int writestep = 1<<10;
    unsigned int numberOfAcceptedSteps;
    for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
        numberOfAcceptedSteps = 0;
        for (unsigned int j = 0; j < m_numberOfParticles; j++) {
            // Call solver method to do a single Monte-Carlo step.
            bool acceptedStep = m_solver->step(timestep, *m_waveFunction, m_particles, j);

            // Truthiness nonsense to count the number of particle movements that were accepted
            numberOfAcceptedSteps += acceptedStep;
        }
        
        // Avoid race conditions
        #pragma omp critical 
        {
        sampler->sample(numberOfAcceptedSteps, this);
        if (omp_get_thread_num() == 0 && i%writestep == 0) {myfile << sampler->getCumulativeEnergy()/sampler->getStepNumber() << std::endl;}
        }
        }
    }

    myfile.close();
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
