#pragma once

#include "materialparameters.h"
#include "geometry.h"
#include "elasticity.h"

namespace phonon_pumping {
    void fit_sound_velocities(const MaterialParameters& material, const Geometry& geometry, const AcousticFrequencyTable& target_freqs, const enhance_mechanism i_mechanism);
}