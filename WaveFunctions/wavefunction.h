#pragma once
#include <memory>
#include <vector>


class WaveFunction {
public:
    // Use the default deconstructor generated by the compiler
    virtual ~WaveFunction() = default;

    int getNumberOfParameters() { return m_numberOfParameters; }
    const std::vector<double>& getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<std::unique_ptr<class Particle>>& particles) = 0;
    virtual double evaluate1D(double x, unsigned int dimension) = 0;
    virtual double reducedF(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k) = 0;
    virtual std::vector<double> quantumForce1D(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k) = 0;
    virtual double wfDerivative(std::vector<std::unique_ptr<class Particle>>& particles) = 0;
    virtual double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) = 0;
    virtual void adjustAlpha(double adjust) = 0;

protected:
    int m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    double m_a = 0;
};

