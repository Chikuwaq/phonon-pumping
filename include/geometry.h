#include <cstdio>
#include <stdexcept>
#include <iostream>

#include "vector.h"
#include "parser.h"
#include "materialparameters.h"

namespace phonon_pumping {

    struct Geometry {
        private:
        static uint16_t expect_uint(const double value);
        static std::tuple<uint16_t, uint16_t> get_index_range(std::map<std::string, double>& parsed_input, const std::string& key);

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

        // label-value pairs of requested output
        std::vector<std::pair<std::string, double>> horizontal_lines;

        void load(const std::filesystem::path& filepath, const MaterialParameters& material);

        void setup_grids();
    };

}