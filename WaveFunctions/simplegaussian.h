#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(double alpha);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double evaluate1D(double x, unsigned int dimension);
    double reducedF(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k);
    std::vector<double> quantumForce1D(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k);
    double wfDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double geta() {return m_a;};
    void adjustAlpha(double d);
};
