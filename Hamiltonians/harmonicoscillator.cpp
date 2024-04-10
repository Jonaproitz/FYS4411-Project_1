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
    std::vector<double> r_2 = particles.at(1)->getPosition();
    double wf = 1;
    for (unsigned int i=0; i<r_2.size(); i++) {
        wf *= waveFunction.evaluate1D(r_2.at(i), i);
    }
    wf *= waveFunction.reducedF(particles, r_2, 1);

    return 2*wf*wf;
}