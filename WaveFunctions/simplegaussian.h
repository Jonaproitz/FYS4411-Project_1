#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(double alpha);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double evaluate1D(double x, unsigned int dimension);
    double reducedF(std::vector<std::unique_ptr<class Particle>>& particles, unsigned int particle_i, unsigned int dimension, double sl);
    double quantumForce1D(std::vector<std::unique_ptr<class Particle>>& particles, unsigned int particle_i, double x, unsigned int dimension, double sl);
    double wfDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    void adjustAlpha(double d);
};
