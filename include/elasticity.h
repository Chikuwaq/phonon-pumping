/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include "materialparameters.h"
#include "geometry.h"

namespace phonon_pumping {

    void calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw);
    
    double eigenvalue_equation(const MaterialParameters& material, const Geometry& geometry, const bool is_TA, const double freq_GHz);
    void find_acoustic_eigenfrequencies(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw);
    void append_eigenfrequencies(const MaterialParameters& material, const Geometry& geometry, const enhance_mechanism i_mechanism, list_of_frequencies& frequencies_to_draw);
}