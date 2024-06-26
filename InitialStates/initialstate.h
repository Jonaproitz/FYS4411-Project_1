#pragma once

#include <memory>
#include <vector>

#include "../particle.h"
#include "Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double timestep,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& randomEngine
            );
