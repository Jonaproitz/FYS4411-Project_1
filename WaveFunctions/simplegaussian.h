#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(double alpha);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double evaluate1D(double x);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    void adjustAlpha(double d);
};
