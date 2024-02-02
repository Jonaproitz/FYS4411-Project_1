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
        std::vector<std::unique_ptr<class Particle>>& particles){
    Random ran;
    double sl = (ran.nextDouble()-0.5)*stepLength;
    double val = ran.nextDouble();
    double x = particles[0]->getPosition()[0] + sl;
    double wav = exp(-(0.5*x*x));
    double wav_old = waveFunction.evaluate(particles);
    double check = wav*wav / (wav_old*wav_old);
    bool a = val<=check;

    if (a)
    {particles[0]->adjustPosition(sl, 0);
}

    return a;
}
