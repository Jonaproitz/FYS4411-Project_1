#include <memory>
#include <cmath>
#include <cassert>
#include <iostream>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    double alpha = m_parameters[0];
    double x = particles[0]->getPosition()[0];
    return exp(-alpha*(x*x));
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
    double alpha = m_parameters[0];
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
