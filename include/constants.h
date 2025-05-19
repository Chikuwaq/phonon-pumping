/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

namespace constants {
    constexpr double PI = 3.141592653589;

    constexpr double MU0 = 4.0 * PI * 1e-7;  // [N A^(-2)] permeability in vacuum
    constexpr double GAMMA = 2.0 * PI * 28.5e9;  // [rad.Hz/T] [An2020] electron gyromagnetic ratio
    constexpr double SCALE_TO_GIGA = 1e-9;
    constexpr double SCALE_TO_NANO = 1e9;
};

namespace phonon_pumping {

    enum enhance_mechanism {
        STRESS_MATCHING_TA_NeumannNeumann,
        STRESS_MATCHING_LA_NeumannNeumann,
        STRESS_MATCHING_TA_DirichletNeumann,
        EIGENMODE_TA,
        EIGENMODE_LA, // Check: this resonance exists, doesn't it?
        NUM_ELEMENTS
    };

};

typedef std::array<std::vector<std::pair<std::string, double>>, phonon_pumping::enhance_mechanism::NUM_ELEMENTS> list_of_frequencies;