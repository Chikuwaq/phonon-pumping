#pragma once

#include "materialparameters.h"
#include "geometry.h"

namespace phonon_pumping {

    void calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw);
    double lhs_Eq38(const MaterialParameters& material, const Geometry& geometry, const double freq_GHz);
    void find_acoustic_eigenfrequencies(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw);

}