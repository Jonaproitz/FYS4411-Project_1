#include<memory>
#include <cassert>
#include <iostream>
#include <vector>
#include <math.h>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(double omega){
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles){
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++) {
        std::vector<double> r = particles[i]->getPosition();
        for (unsigned int di = 0; di < r.size(); di++){
            double x = r.at(di);
            if (di == 2) {r2 += waveFunction.getParameters().at(1)*waveFunction.getParameters().at(1)*x*x;}
            else {r2 += x*x;}
        }
    }
    return 0.5*(-waveFunction.computeDoubleDerivative(particles) + r2);
}

double HarmonicOscillator::computeOnebodyDensity(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles) {
    double a = waveFunction.geta();
    std::vector<double> r_1 = particles.at(0)->getPosition();
    double alpha = waveFunction.getParameters().at(0);
    if (a != 0) {
        std::vector<double> r_2 = particles.at(1)->getPosition();
        double u2 = 0, r2 = 0;  
        for (unsigned int i=0; i<r_1.size(); i++) {
            u2 += (r_2.at(i) - r_1.at(i))*(r_2.at(i) - r_1.at(i));
            r2 += r_1.at(i)*r_1.at(i);
        }
        
        double u = sqrt(u2);
        double r = sqrt(r2);
        
        return 1/(2*M_PI) / (alpha * r) * exp(-4*alpha*r2) * u * sinh(4*alpha*u*r)*exp(-2*alpha*u)*(1 - a/u)*(1 - a/u);
    }
    else {
        double r2 = 0;
        for (double x:r_1) {
            r2 += x*x;
        }
        return sqrt(8 * alpha / M_PI) * exp(-2 * alpha * r2);
    }

}