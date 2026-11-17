/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string_view>
#include "materialparameters.h"
#include "geometry.h"

namespace phonon_pumping {

    const std::string spectrum_datafile_name = "spectrum.dat";
    const std::string gnuplot_file_name = "plot_spectrum.plt";

    // Approximation of the perceptually uniform, colorblind-friendly
    // "cividis" colormap, sampled at 6 points and interpolated linearly
    // by gnuplot.
    constexpr std::string_view gnuplot_palette_cividis() {
        return "set palette defined ("
               "0 '#00204c', "
               "1 '#31446b', "
               "2 '#666970', "
               "3 '#958f78', "
               "4 '#cbba69', "
               "5 '#ffe945')\n";
    }

    void write_gnuplot_script(const MaterialParameters& material, const Geometry& geometry, const list_of_frequencies& frequencies_to_draw);
    void write_gnuplot_command_draw_frequency_line(const MaterialParameters& material, const Geometry& geometry, const std::string& color_name, const std::vector<std::pair<std::string, double>>& frequencies_to_draw, const bool outside, std::ofstream& file);
    double screen_position_x(const Geometry& geometry, const double position_with_respect_to_screen);

    void write_heatmap_xy(const std::vector<std::vector<double>>& data, const std::vector<double>& x, const std::vector<double>& y, const std::string& x_label, const std::string& y_label, const std::string& cb_label, const std::string& data_file, const std::string& script_file, const std::string& image_file);

}