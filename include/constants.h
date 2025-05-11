#pragma once

namespace constants {
    constexpr double PI = 3.141592653589;

    constexpr double MU0 = 4.0 * PI * 1e-7;  // [N A^(-2)] permeability in vacuum
    constexpr double GAMMA = 2.0 * PI * 28.5e9;  // [rad.Hz/T] [An2020] electron gyromagnetic ratio
    constexpr double SCALE_TO_GIGA = 1e-9;
    constexpr double SCALE_TO_NANO = 1e9;
};
