#include <cstdio>
#include <stdexcept>
#include <iostream>

#include "vector.h"
#include "parser.h"

namespace phonon_pumping {

    struct Geometry {
        private:

        public:
        // geometry
        double mag_thickness_nm;
        double nonmag_thickness_nm;
        double theta;

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

        void read(const std::filesystem::path& filepath) {
            std::cout << "\nReading geometry file..." << std::endl;

            auto result = parse_file(filepath);

            if (result.find("theta") == result.end()) {
                throw std::runtime_error("Magnetization angle is missing in the input file");
            }
            if (result["theta"] != 0.0) {
                throw std::runtime_error("Nonzero magnetization angle is not yet implemented");
            }
            for (const auto& [key, value] : result) {
                std::cout << key << " = " << value << std::endl;
            }

            mag_thickness_nm = result["mag_thickness_nm"];
            nonmag_thickness_nm = result["nonmag_thickness_nm"];
            theta = result["theta"];

            h_static = result["h_static"];
            hpara2 = result["hpara2"];
            hperp2 = result["hperp2"];

            min_B_Tesla = result["min_B_Tesla"];
            max_B_Tesla = result["max_B_Tesla"];
            n_B = result["n_B"];
            min_freq_GHz = result["min_freq_GHz"];
            max_freq_GHz = result["max_freq_GHz"];
            n_freq = result["n_freq"];
        }

        void setup_grids() {
            B_Tesla.linspace(min_B_Tesla, max_B_Tesla, n_B);
            freq_GHz.linspace(min_freq_GHz, max_freq_GHz, n_freq);
        }
    };

}