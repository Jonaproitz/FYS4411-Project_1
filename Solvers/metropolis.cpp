#include <memory>
#include <vector>

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
        unsigned int particle_i,
        unsigned int dimension){
    double D = 0.5;
    double x = particles.at(particle_i)->getPosition().at(dimension);
    double wf_old = waveFunction.evaluate1D(x, dimension)*waveFunction.reducedF(particles, particle_i, dimension, 0);
    double qf_old = waveFunction.quantumForce1D(particles, particle_i, x, dimension, 0);

    double sl = m_rng->nextGaussian(0, 1)*sqrt(timestep) + qf_old*timestep*D;
    double x_new = x+sl;
    double wf_new = waveFunction.evaluate1D(x_new, dimension)*waveFunction.reducedF(particles, particle_i, dimension, sl);
    double qf_new = waveFunction.quantumForce1D(particles, particle_i, x_new, dimension, sl);

    double Greensfunction = 0.5*(qf_old+qf_new)*(D*timestep*0.5*(qf_old-qf_new) - x_new+x);
    Greensfunction = exp(Greensfunction);

    double check = Greensfunction * wf_new*wf_new/(wf_old*wf_old);
    bool a = m_rng->nextDouble()<=check;

    if (a == true)
    {particles[particle_i]->adjustPosition(sl, dimension);
}
    
    return a;
}
