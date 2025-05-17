#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <array>

#include "vector.h"
#include "parser.h"
#include "materialparameters.h"

namespace phonon_pumping {

    enum Enhance_mechanism {
        STRESS_MATCHING_TA_NeumannNeumann,
        STRESS_MATCHING_LA_NeumannNeumann,
        STRESS_MATCHING_TA_DirichletNeumann,
        NUM_ELEMENTS
    };

    struct Geometry {
        private:
        static uint16_t expect_uint(const double value);
        
        public:
        // geometry
        double mag_thickness_nm;
        double nonmag_thickness_nm;
        double magnetization_angle_radian;

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

        // label-value pairs of frequency line output
        std::array<std::vector<std::pair<std::string, double>>, Enhance_mechanism::NUM_ELEMENTS> horizontal_lines;

        void load(const std::filesystem::path& filepath, const MaterialParameters& material);

        void setup_grids();
    };

}