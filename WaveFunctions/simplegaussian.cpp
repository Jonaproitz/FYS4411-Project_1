#include <memory>
#include <cmath>
#include <cassert>
#include <iostream>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

// Define parameters of the wavefunction
SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    // Define the wavefunction
    double alpha = m_parameters.back();
    double E = 1;
    double p = particles.size();
    double d = particles[0]->getPosition().size();
    for (unsigned int i = 0; i < p; i++){
        for (unsigned int j = 0; j < d; j++) {
            double x = particles.at(i)->getPosition().at(j);
            E *= exp(-alpha*(x*x));
        }
    }
    return E;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
    // Compute the local kinetic part of the local derivative
    double alpha = m_parameters.at(0);
    double E = 0.0;
    for (unsigned int j=0; j < particles.size(); j++) {
        std::vector<double> r = particles[j]->getPosition();
        double r2 = 0;
        unsigned int d = r.size();
        for (unsigned int i=0; i<d; i++){
            r2 = r2 + r[i]*r[i];
        }
        E += 2*alpha*(2*alpha*r2 - d);
    }
    return E;
}

void SimpleGaussian::adjustAlpha(double adjust) {
    // Adjust alpha, and store the new value
    m_numberOfParameters += 1;
    m_parameters.at(0) += adjust;
}