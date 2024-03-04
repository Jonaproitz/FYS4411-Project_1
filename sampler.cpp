#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

#include <iomanip>

using std::cout;
using std::endl;


Sampler::Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double timestep,
        unsigned int numberOfMetropolisSteps)
{
    m_stepNumber = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    m_cumulativeEnergy = 0;
    m_cumulativeEnergy2 = 0;
    m_DerivativePsiE = 0;
    m_DeltaPsi = 0;
    m_EnergyDer = 0;
    m_stepLength = timestep;
    m_numberOfAcceptedSteps = 0;
}


void Sampler::sample(unsigned int acceptedStep, System* system) {
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2  += localEnergy*localEnergy;

    auto WFder = system->wfDerivative();
    m_DeltaPsi += WFder;
    m_DerivativePsiE += WFder*localEnergy;

    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;
}

void Sampler::printOutputToTerminal(System& system) {
    auto pa = system.getWaveFunctionParameters();
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(m_numberOfMetropolisSteps) << endl;
    cout << " Step length used : " << m_stepLength << endl;
    cout << " Ratio of accepted steps: " << ((double) m_numberOfAcceptedSteps) / ((double) m_numberOfMetropolisSteps*m_numberOfParticles*m_numberOfDimensions*m_energy.size()) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    cout << endl;
    cout << "  -- Results --  " << endl;
    const char separator    = ' ';
    const int nameWidth     = 16;
    cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Alpha values";
    cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Energies";
    cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Variance";
    cout << endl;
    for (unsigned int i=0; i<m_energy.size(); i++){
        cout << std::left << std::setw(nameWidth) << std::setfill(separator) << m_alphaValues.at(i);
        cout << std::left << std::setw(nameWidth) << std::setfill(separator) << m_energy.at(i);
        cout << std::left << std::setw(nameWidth) << std::setfill(separator) << m_variance.at(i);
        cout << endl;
    }
}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities.
    m_energy.push_back(m_cumulativeEnergy / m_numberOfMetropolisSteps);
    double m_energy2 = m_cumulativeEnergy2 / m_numberOfMetropolisSteps;
    m_variance.push_back(m_energy2 - m_energy.back()*m_energy.back());

    m_EnergyDer = 2*(m_DerivativePsiE-m_DeltaPsi*m_energy.back())/m_numberOfMetropolisSteps;

    // Reset cumulative properties for next run
    m_cumulativeEnergy = 0;
    m_cumulativeEnergy2 = 0;

    m_DerivativePsiE = 0;
    m_DeltaPsi = 0;
}
