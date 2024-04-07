#include <memory>
#include <iostream>
#include <cassert>

#include "initialstate.h"
#include "../particle.h"
#include "../Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double timestep,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        // Give every particle a position in space chosen from a gaussian distribution in each dimension
        for (unsigned int j=0; j < numberOfDimensions; j++) {
            double num = rng.nextGaussian(0, 1);
            position.push_back(num*sqrt(timestep));
        }

        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}
