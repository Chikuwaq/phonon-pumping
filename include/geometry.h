/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <stdexcept>
#include <iostream>
#include <array>

#include "vector.h"
#include "parser.h"
#include "materialparameters.h"

namespace phonon_pumping {

    struct Bilayer {
        double mag_thickness_nm;
        double nonmag_thickness_nm;

        public:
        Bilayer(const double d, const double L) 
        : mag_thickness_nm(d), nonmag_thickness_nm(L) {};

    };

    struct Geometry {
        private:
        static uint16_t expect_uint(const double value);
        
        public:
        // geometry
        Bilayer bilayer{0.0, 0.0};
        double magnetization_angle_radian;

        double h_static;
        // ac microwave intensity
        double hpara2;
        double hperp2;

        // sweep range
        double min_B_Tesla;
        double max_B_Tesla;
        size_t n_B;
        double min_freq_GHz;
        double max_freq_GHz;
        size_t n_freq;

        // calculated Vectors
        Vector<double> B_Tesla;
        Vector<double> freq_GHz;

        std::array<bool, enhance_mechanism::NUM_ELEMENTS> should_draw_frequencies;

        public:
        void load(const std::filesystem::path& filepath);
        void setup_grids();
    };

}