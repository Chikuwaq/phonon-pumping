#pragma once

#include "materialparameters.h"
#include "geometry.h"

namespace phonon_pumping {

    void calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw);

}