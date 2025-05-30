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

    const std::string spectrum_datafile_name = "spectrum.dat";
    const std::string gnuplot_file_name = "plot_spectrum.plt";

    void output_spectrum(const Vector<double>& x_values, const Vector<double>& y_values, const Matrix2D<double>& z_values);

    void write_gnuplot_script(const MaterialParameters& material, const Geometry& geometry, const list_of_frequencies& frequencies_to_draw);
    void write_gnuplot_command_draw_frequency_line(const MaterialParameters& material, const Geometry& geometry, const std::string& color_name, const std::vector<std::pair<std::string, double>>& frequencies_to_draw, const bool outside, std::ofstream& file);
    double screen_position_x(const Geometry& geometry, const double position_with_respect_to_screen);

}