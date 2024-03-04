#include<memory>
#include <cassert>
#include <iostream>
#include <vector>

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
        for (double x:r){
            r2 += x*x;
        }
    }
    return 0.5*(-waveFunction.computeDoubleDerivative(particles) + m_omega*m_omega*r2);
}
