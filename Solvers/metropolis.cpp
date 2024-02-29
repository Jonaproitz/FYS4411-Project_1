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
        double stepLength,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles,
        unsigned int particle_i,
        unsigned int dimension){
    double alpha = waveFunction.getParameters().at(0);
    Random ran;
    double sl = (ran.nextDouble()-0.5)*stepLength;
    double x = particles.at(particle_i)->getPosition().at(dimension);
    double wf_old = waveFunction.evaluate1D(x);
    double wf_new = waveFunction.evaluate1D(x+sl);
    double check = wf_new*wf_new/(wf_old*wf_old);
    bool a = ran.nextDouble()<=check;

    if (a == true)
    {particles[particle_i]->adjustPosition(sl, dimension);
}
    
    return a;
}
