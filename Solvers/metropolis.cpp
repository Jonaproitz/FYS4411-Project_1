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
    double alpha = waveFunction.getParameters()[0];
    Random ran;
    double sl = (ran.nextDouble()-0.5)*stepLength;
    double x = particles[particle_i]->getPosition()[dimension];
    double check = (2*alpha*(x+sl)*(x+sl) - 1)*(2*alpha*(x+sl)*(x+sl) - 1)/((2*alpha*x*x - 1)*(2*alpha*x*x - 1));
    bool a = ran.nextDouble()<=check;

    if (a)
    {particles[particle_i]->adjustPosition(sl, dimension);
}
    
    return a;
}
