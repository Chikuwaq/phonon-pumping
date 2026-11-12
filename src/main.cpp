/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <filesystem>

#include "constants.h"
#include "matrix.h"
#include "materialparameters.h"
#include "geometry.h"
#include "elasticity.h"
#include "magnetoelasticity.h"
#include "gnuplot.h"
#include "task.h"

using namespace phonon_pumping;

constexpr std::string_view VERSION_NO = "1.0.1";

void print_usage() {
    std::cout << "Usage: ./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file -t path/to/task_file" << std::endl;
}

void plot_FMR_spectrum(const MaterialParameters& material, const Geometry& geometry) {
    // label-value pairs of frequency line output
    list_of_frequencies frequencies_to_draw;
    calculate_Gilbert_damping_enhancing_freq(material, geometry, frequencies_to_draw);
    find_acoustic_eigenfrequencies(material, geometry, frequencies_to_draw);

    // --- compute FMR absorption spectrum
    Matrix2D<double> FMR_power_absorption; 
    calculate_spectrum(material, geometry, FMR_power_absorption);
    output_spectrum(geometry.B_Tesla, geometry.freq_GHz, FMR_power_absorption);

    write_gnuplot_script(material, geometry, frequencies_to_draw);
    std::cout << std::endl;
    std::cout << "Run the gnuplot script " << gnuplot_file_name << " to generate figures." << std::endl;
}

void fit_sound_velocities(const MaterialParameters& material, const Geometry& geometry, const std::vector<double>& target_freqs) {
    // copy data not to damage the original
    MaterialParameters materialTemp = material;

    constexpr double n_sound_velocities = 100.0;
    Vector<double> sound_velocities_magnet; 
    sound_velocities_magnet.linspace(material.min_cl_ms, material.max_cl_ms, n_sound_velocities);
    Vector<double> sound_velocities_nonmagnet = sound_velocities_magnet;

    Matrix2D<double> deviations;
    deviations.zeros(n_sound_velocities, n_sound_velocities);

    for (uint16_t i_sound_velocity_magnet = 0; i_sound_velocity_magnet < n_sound_velocities; i_sound_velocity_magnet++) {
        materialTemp.mag_cl_ms = sound_velocities_magnet(i_sound_velocity_magnet);

        for (uint16_t i_sound_velocity_nonmagnet = 0; i_sound_velocity_nonmagnet < n_sound_velocities; i_sound_velocity_nonmagnet++) {
            materialTemp.nonmag_cl_ms = sound_velocities_nonmagnet(i_sound_velocity_nonmagnet);

            list_of_frequencies eigenfrequencies;
            find_acoustic_eigenfrequencies(materialTemp, geometry, eigenfrequencies);

            std::vector<double> LA_frequencies = extract_LA_frequencies(eigenfrequencies);
            const uint16_t n_freqs = std::min(target_freqs.size(), LA_frequencies.size());
            double sum = 0.0;
            for (uint16_t i_freq = 0; i_freq < n_freqs; i_freq++) {
                const double deviation = target_freqs[i_freq] - LA_frequencies[i_freq];
                sum += deviation * deviation;
            }
            std::cout << "sqrt(Sum of squared frequency deviation) = " << std::sqrt(sum) << std::endl;
            deviations(i_sound_velocity_magnet, i_sound_velocity_nonmagnet) = std::sqrt(sum);
        }
    }
    auto [irow, icol, min_value] = deviations.min();
    std::cout << "Fitting results: " << std::endl;
    std::cout << "\tSound velocity of magnet = " << sound_velocities_magnet(irow) << std::endl;
    std::cout << "\tSound velocity of nonmagnet = " << sound_velocities_nonmagnet(icol) << std::endl;
    std::cout << "\tsqrt(Sum of squared frequency deviation) = " << min_value << std::endl;
    

}


int main(int argc, char *argv[]) {

    if (argc < 5) {
        print_usage();
        return -1;
    }

    std::filesystem::path path_material_parameter;
    std::filesystem::path path_geometry;
    std::filesystem::path path_task;

    for(int count = 1; count < argc; count++ ) {
        std::string s = argv[count];
        if (s == "-m") {
            count++;
            path_material_parameter = argv[count];
        }
        else if (s == "-g") {
            count++;
            path_geometry = argv[count];
        }
        else if (s == "-t") {
            count++;
            path_task = argv[count];
        }
    }

    if (path_material_parameter.empty()) {
        std::cout << "Path to material parameter file is empty!" << std::endl;
        print_usage();
        return -1;
    }
    if (path_geometry.empty()) {
        std::cout << "Path to geometry file is empty!" << std::endl;
        print_usage();
        return -1;
    }

    try {
        MaterialParameters material;
        material.set(path_material_parameter);

        Geometry geometry;
        geometry.load(path_geometry);
        geometry.setup_grids();

        std::vector<task::TaskID> tasks;
        if (path_task.empty()) {
            std::cout << "Path to simulation control file is empty. Fall back to default (plot FMR spectrum)..." << std::endl;
            tasks.push_back(task::TaskID::PLOT_FMR_SPECTRUM);
        }
        else {
            tasks = task::load(path_task);
        }

        for (auto& task : tasks) {
            if (task == task::TaskID::PLOT_FMR_SPECTRUM) {
                plot_FMR_spectrum(material, geometry);
            }
            else if (task == task::TaskID::FIT_SOUND_VELOCITIES) {
                std::vector<double> target_freqs = task::extract_freqs(path_task);
                fit_sound_velocities(material, geometry, target_freqs);
            }
        }

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "DONE. " << std::endl;
    std::cout << "\n-----------------------------------------------------------------" << std::endl;
    std::cout << "Phonon-pumping plotter " << VERSION_NO << " © 2025 Takuma Sato" << std::endl;
    std::cout << "This program is distributed under the GNU General Public License \nversion 3.0, or any later version." << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    
    return 0;
}