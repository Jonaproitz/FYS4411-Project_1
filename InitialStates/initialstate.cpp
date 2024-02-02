#include <memory>
#include <iostream>
#include <cassert>

#include "initialstate.h"
#include "../particle.h"
#include "../Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double stepLength,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (unsigned int j=0; j < numberOfDimensions; j++) {
            Random ran;
            double num = ran.nextDouble()-0.5;
            position.push_back(num);
        }

        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}
