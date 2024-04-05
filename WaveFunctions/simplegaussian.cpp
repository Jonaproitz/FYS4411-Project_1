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
    m_a = 0.0043;
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    // Define the wavefunction in nD for n particles from the definition of 1D
    // This function isnt really used
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

double SimpleGaussian::reducedF(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k) {
    // Define the interaction term for a single particle, interacting with all other particles
    unsigned int p = particles.size();
    unsigned int d = rk.size();
    double f = 1;
    for (unsigned int j=0; j<p; j++) {
        if (j != particle_k && m_a != 0) {
            std::vector<double> rj = particles.at(j)->getPosition();
            double rkl, rkl2 = 0;
            for (unsigned int di=0; di<d; di++) {
                rkl2 = (rj.at(di) - rk.at(di))*(rj.at(di) - rk.at(di));
            }
            rkl = sqrt(rkl2);
            if (rkl <= m_a) {f = 1e-10; break;}
            else {f *= 1 - m_a/rkl;}
        }
    }
    return f;
}

std::vector<double> SimpleGaussian::quantumForce1D(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double> rk, unsigned int particle_k) {
    double beta;

    std::vector<double> qf;
    unsigned int dimension = rk.size();
    for (unsigned int dim = 0; dim < dimension; dim++) {
        qf.push_back(0.0);
    }
    for (unsigned int j=0; j<particles.size(); j++) {
        if (j != particle_k && m_a != 0) {
            std::vector<double> rj = particles.at(j)->getPosition();
            // Find the length between particles
            double rkj2 = 0.0;
            for (unsigned int di=0; di<rk.size(); di++) {
                double x = rk.at(di) - rj.at(di);
                rkj2 += x*x;
            }
            double rkj = sqrt(rkj2);
            if (rkj <= m_a) {
                for (unsigned int i = 0; i<dimension; i++) {
                    qf.at(i) = (rk.at(i) - rj.at(i))*1e10;
                }
                break;
            }
            for (unsigned int dim=0; dim<dimension; dim++) {
                // Find the interaction between particles
                qf.at(dim) += (rk.at(dim) - rj.at(dim)) * m_a/(sqrt(rkj2) * (1- m_a/sqrt(rkj2)));
            }
        }
    }


    for (unsigned int i=0; i < dimension; i++) {
        if (i == 2) {beta = m_parameters.at(1);}
        else {beta = 1.0;}
        // Store the final expression for the quantum force
        qf.at(i) = -4*m_parameters.at(0)*beta*rk.at(i) + 2*qf.at(i);
    }
    return qf;
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
            if (j != k && m_a != 0) {
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
                E += -4*alpha * temp*m_a/(rkj*rkj*rkj*(1-m_a/rkj));
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
                        E += temp/(rkj*rki) * m_a/(rkj2*(1-m_a/rkj)) * m_a/(rki2*(1 - m_a/rki));
                    }
                }
                // Fourth term
                E += 2*m_a/(rkj*rkj2*(1-m_a/rkj));
                E += m_a/(rkj2*rkj*(1 - m_a/rkj)) * (m_a/(rkj*(1-m_a/rkj)) - 1);
            }
        }
    }
    return E;
}

void SimpleGaussian::adjustAlpha(double adjust) {
    // Adjust alpha value
    m_parameters.at(0) += adjust;
}