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
    /*
    std::vector<double> r = particles[0]->getPosition();
    for (unsigned int i = 0; i<r.size(); i++){
        double x = r[i];
        double x_new = x + sl;
        double check = exp(alpha*(x-x_new));
        bool a = val<=check;
        if (a)
            {particles[0]->adjustPosition(sl, i);}
    }
    */
    //double x = particles[0]->getPosition()[dimension];
    //double wav = exp(-(0.5*x*x));
    //double wav_old = waveFunction.evaluate(particles);
    //double check = wav*wav / (wav_old*wav_old);
    double x = particles[particle_i]->getPosition()[dimension];
    double check = (2*alpha*(2*alpha*(x+sl)*(x+sl) - 1))*(2*alpha*(2*alpha*(x+sl)*(x+sl) - 1))/(2*alpha*(2*alpha*x*x - 1)*2*alpha*(2*alpha*x*x - 1));//exp(2*alpha*sl);
    bool a = ran.nextDouble()<=check;

    if (a)
    {particles[particle_i]->adjustPosition(sl, dimension);
}
    
    return a;
}
