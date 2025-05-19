/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include "matrix.h"
#include "materialparameters.h"
#include "geometry.h"


namespace phonon_pumping {
        
    void calculate_spectrum(const MaterialParameters& material, const Geometry& geometry, Matrix2D<double>& FMR_power_absorption);
    
}
