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
#include "fitting.h"
#include "csv.h"

using namespace phonon_pumping;

constexpr std::string_view VERSION_NO = "1.2.0";

void print_usage() {
    std::cout << "Usage: ./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file -t path/to/task_file -c path/to/frequency_csv_file" << std::endl;
}

void plot_FMR_spectrum(const MaterialParameters& material, const Geometry& geometry) {
    // label-value pairs of frequency line output
    list_of_frequencies frequencies_to_draw;
    calculate_Gilbert_damping_enhancing_freq(material, geometry, frequencies_to_draw);
    find_acoustic_eigenfrequencies(material, geometry, geometry.should_draw_frequencies[enhance_mechanism::EIGENMODE_TA], geometry.should_draw_frequencies[enhance_mechanism::EIGENMODE_LA], frequencies_to_draw);

    // --- compute FMR absorption spectrum
    Matrix2D<double> FMR_power_absorption; 
    calculate_spectrum(material, geometry, FMR_power_absorption);
    FMR_power_absorption.output_spectrum(geometry.B_Tesla, geometry.freq_GHz, spectrum_datafile_name);

    write_gnuplot_script(material, geometry, frequencies_to_draw);
    std::cout << std::endl;
    std::cout << "Run the gnuplot script " << gnuplot_file_name << " to generate figures." << std::endl;
}


int main(int argc, char *argv[]) {

    if (argc < 9) {
        print_usage();
        return EXIT_FAILURE;
    }

    std::filesystem::path path_material_parameter, path_geometry, path_task, path_csv;

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
        else if (s == "-c") {
            count++;
            path_csv = argv[count];
        }
    }

    if (path_material_parameter.empty()) {
        std::cout << "Path to material parameter file is empty!" << std::endl;
        print_usage();
        return EXIT_FAILURE;
    }
    if (path_geometry.empty()) {
        std::cout << "Path to geometry file is empty!" << std::endl;
        print_usage();
        return EXIT_FAILURE;
    }

    try {
        MaterialParameters material;
        material.set(path_material_parameter);

        Geometry geometry;
        geometry.load(path_geometry);
        geometry.setup_grids();

        std::vector<task::TaskID> tasks;
        AcousticFrequencyTable target_freqs_table;

        if (path_task.empty()) {
            std::cout << "Path to simulation control file is empty. Fall back to default (plot FMR spectrum)..." << std::endl;
            tasks.push_back(task::TaskID::PLOT_FMR_SPECTRUM);
        }
        else {
            tasks = task::load(path_task);
            if (!path_csv.empty()) {
                Matrix2D<std::string> target_freq_table;
                const bool success = read_csv(path_csv, target_freq_table) && target_freqs_table.load(target_freq_table);
                if (!success) {
                    std::cout << "Could not read CSV file" << std::endl;
                    return EXIT_FAILURE;
                }
                if (target_freqs_table.has_mode(enhance_mechanism::EIGENMODE_TA)) tasks.push_back(task::TaskID::FIT_TA_SOUND_VELOCITIES);
                if (target_freqs_table.has_mode(enhance_mechanism::EIGENMODE_LA)) tasks.push_back(task::TaskID::FIT_LA_SOUND_VELOCITIES);
            }
        }

        for (auto& task : tasks) {
            if (task == task::TaskID::PLOT_FMR_SPECTRUM) {
                plot_FMR_spectrum(material, geometry);
            }
            else if ((task == task::TaskID::FIT_TA_SOUND_VELOCITIES) || (task == task::TaskID::FIT_LA_SOUND_VELOCITIES)) {
                enhance_mechanism i_mechanism = (task == task::TaskID::FIT_TA_SOUND_VELOCITIES) ? enhance_mechanism::EIGENMODE_TA : enhance_mechanism::EIGENMODE_LA;

                fit_sound_velocities(material, geometry, target_freqs_table, i_mechanism);
            }
            else {
                std::cout << "Task file contains unknown task" << std::endl;
                return EXIT_FAILURE;
            }
        }

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "DONE. " << std::endl;
    std::cout << "\n-----------------------------------------------------------------" << std::endl;
    std::cout << "Phonon-pumping plotter " << VERSION_NO << " © 2025 Takuma Sato" << std::endl;
    std::cout << "This program is distributed under the GNU General Public License \nversion 3.0, or any later version." << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    
    return EXIT_SUCCESS;
}