#include <iostream>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <cmath>
#include <array>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#include "constants.h"
#include "matrix.h"
#include "materialparameters.h"
#include "geometry.h"
#include "elasticity.h"
#include "magnetoelasticity.h"
#include "gnuplot.h"

using namespace phonon_pumping;


void print_usage() {
    std::cout << "Usage: ./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file" << std::endl;
}


int main(int argc, char *argv[]) {

    if (argc < 5) {
        print_usage();
        return -1;
    }

    std::filesystem::path path_material_parameter;
    std::filesystem::path path_geometry;

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

        // label-value pairs of frequency line output
        list_of_frequencies frequencies_to_draw;
        calculate_Gilbert_damping_enhancing_freq(material, geometry, frequencies_to_draw);
        find_acoustic_eigenfrequencies(material, geometry, frequencies_to_draw);

        // --- compute FMR absorption spectrum
        Matrix2D<double> FMR_power_absorption; 
        calculate_spectrum(material, geometry, FMR_power_absorption);
        output_spectrum(geometry.B_Tesla, geometry.freq_GHz, FMR_power_absorption);

        write_gnuplot_script(material, geometry, frequencies_to_draw);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "\nPhononPumpingPlotter DONE. Run the gnuplot script " << gnuplot_file_name << " to generate figures." << std::endl;
    
    return 0;
}