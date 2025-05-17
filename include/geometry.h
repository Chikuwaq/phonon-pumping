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

        void load(const std::filesystem::path& filepath);
        void setup_grids();
    };

}