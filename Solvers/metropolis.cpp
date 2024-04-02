#include <memory>
#include <vector>
#include <iostream>

#include "metropolis.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng)){
}


bool Metropolis::step(
        double timestep,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles,
        unsigned int particle_i){
    double D = 0.5;
    std::vector<double> r = particles.at(particle_i)->getPosition();
    std::vector<double> r_new, sl, qf_old, qf_new;
    double f = waveFunction.reducedF(particles, r, particle_i);
    qf_old = waveFunction.quantumForce1D(particles, r, particle_i);
    double wf_old = 1.0, wf_new = 1.0;
    for (unsigned int i=0; i < r.size(); i++) {
        wf_old *= waveFunction.evaluate1D(r.at(i), i)*f;
        sl.push_back(m_rng->nextGaussian(0, 1)*sqrt(timestep) + qf_old.at(i)*timestep*D);
        r_new.push_back(r.at(i) + sl.at(i));
    }
    qf_new = waveFunction.quantumForce1D(particles, r_new, particle_i);
    double f_new = waveFunction.reducedF(particles, r_new, particle_i);
    double Greensfunction = 0.0;
    for (unsigned int i=0; i<r.size(); i++) {
        wf_new *= waveFunction.evaluate1D(r_new.at(i), i)*f_new;
        Greensfunction += 0.5*(qf_old.at(i)+qf_new.at(i))*(D*timestep*0.5*(qf_old.at(i)-qf_new.at(i)) - r_new.at(i)+r.at(i));
    }
    Greensfunction = exp(Greensfunction);
    
    double check = Greensfunction * wf_new*wf_new/(wf_old*wf_old);
    bool a = m_rng->nextDouble()<=check;

    if (a == true) {
        for (unsigned int i = 0; i < r.size(); i++) {
            particles[particle_i]->adjustPosition(sl.at(i), i);
        }
    }
    
    return a;
}
