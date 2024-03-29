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
    double beta = 2.82843;
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    // Define the wavefunction in nD for n particles from the definition of 1D
    double E = 1;
    double p = particles.size();
    double d = particles.at(0)->getPosition().size();
    for (unsigned int i = 0; i < p; i++){
        for (unsigned int j = 0; j < d; j++) {
            double x = particles.at(i)->getPosition().at(j);
            E *= evaluate1D(x, j);
        }
    }
    return E;
}

double SimpleGaussian::evaluate1D(double x, unsigned int dimension) {
    // Define the wavefunction for 1 dimension
    // This is used for optimization in the metroplis algo
    if (dimension == 2) {return exp(-m_parameters.at(0)*(m_parameters.at(1)*x*x));}
    return exp(-m_parameters.at(0)*(x*x));
}

double SimpleGaussian::reducedF(std::vector<std::unique_ptr<class Particle>>& particles, unsigned int particle_i, unsigned int dimension, double sl) {
    double a = 0.0043;
    std::vector<double> rk = particles.at(particle_i)->getPosition();
    rk.at(dimension) += sl;
    unsigned int p = particles.size();
    unsigned int d = rk.size();
    double f = 1;
    for (unsigned int j=0; j<p; j++) {
        if (j != particle_i) {
            std::vector<double> rj = particles.at(j)->getPosition();
            double rkl, rkl2 = 0;
            for (unsigned int di=0; di<d; di++) {
                rkl2 = (rj.at(di) - rk.at(di))*(rj.at(di) - rk.at(di));
            }
            rkl = sqrt(rkl2);
            if (rkl <= a) {f = 0; break;}
            else {f *= 1 - a/rkl;}
        }
    }
    return f;
}

double SimpleGaussian::quantumForce1D(std::vector<std::unique_ptr<class Particle>>& particles, unsigned int particle_i, double x, unsigned int dimension, double sl) {
    double beta, a = 0.0043;
    if (dimension == 2) {beta = m_parameters.at(1);}
    else {beta = 1.0;}

    std::vector<double> rk = particles.at(particle_i)->getPosition();
    rk.at(dimension) += sl;
    double temp = 0.0;
    for (unsigned int j=0; j<particles.size(); j++) {
        if (j != particle_i) {
            std::vector<double> rj = particles.at(j)->getPosition();
            double rkj2 = 0.0;
            for (unsigned int di=0; di<rk.size(); di++) {
                double x = rk.at(di) - rj.at(di);
                rkj2 += x*x;
            }
            temp += (rk.at(dimension) - rj.at(dimension)) * a/(sqrt(rkj2) * (1- a/sqrt(rkj2)));
        }
    }

    return -4*m_parameters.at(0)*beta*x + 2*temp;
}

double SimpleGaussian::wfDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
    double WFder = 0;
    double p = particles.size();
    double d = particles.at(0)->getPosition().size();
    for (unsigned int i = 0; i < p; i++){
        std::vector<double> ri = particles.at(i)->getPosition();
        for (unsigned int j = 0; j < d; j++) {
            double x = ri.at(j);
            if (j == 2) {WFder += -m_parameters.at(1)*x*x;}
            else {WFder += -x*x;}
        }
    }
    return WFder;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
    // Compute the local kinetic part of the local derivative
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);
    double a = 0.0043;
    double E = 0.0;

    unsigned int p = particles.size();
    unsigned int d = particles.at(0)->getPosition().size();

    for (unsigned int k=0; k<p; k++) {
        std::vector<double> rk = particles.at(k)->getPosition();
        

        // First term
        double r2 = 0.0;
        for (unsigned int di=0; di<d; di++) {
            double x = rk.at(di);
            if (di == 2) {r2 += beta*x*x;}
            else {r2 += x*x;}
        }
        E += 2*alpha*(2*alpha*r2 - d);

        for (unsigned int j=0; j<p; j++) {
            if (j != k) {
                std::vector<double> rj = particles.at(j)->getPosition();
                double temp = 0.0, rkj2 = 0.0;
                for (unsigned int di=0; di<d; di++) {
                    double xk = rk.at(di), xj = rj.at(di);
                    double xkj = xk - xj;
                    rkj2 += xkj*xkj;
                    if (di == 2) {temp += beta * xk*xkj;}
                    else {temp += xk*xkj;}
                }
                double rkj = sqrt(rkj2);
                // Second term
                E += -4*alpha * temp*a/(rkj*rkj*rkj*(1-a/rkj));
                for (unsigned int i=0; i<p; i++) {
                    if (i != k) {
                        double temp = 0.0, rki2 = 0.0;
                        std::vector<double> ri = particles.at(i)->getPosition();
                        for (unsigned int di=0; di<d; di++) {
                            double xi = ri.at(di);
                            double xk = rk.at(di);
                            temp += (xk - rj.at(di))*(xk - xi);
                            rki2 += (xk-xi)*(xk-xi);
                        }
                        double rki = sqrt(rki2);
                        // Third term
                        E += temp/(rkj*rki) * a/(rkj2*(1-a/rkj)) * a/(rki2*(1 - a/rki));
                    }
                }
                // Fourth term
                E += 2*a/(rkj*rkj2*(1-a/rkj));
                E += a/(rkj2*rkj*(1 - a/rkj)) * (a/(rkj*(1-a/rkj)) - 1);
            }
        }
    }
    return E;
}

void SimpleGaussian::adjustAlpha(double adjust) {
    // Adjust alpha value
    m_parameters.at(0) += adjust;
}