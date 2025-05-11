#include <cstdio>
#include <stdexcept>

#include "vector.h"

namespace phonon_pumping {

    struct Inputs {
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

        void read() {
            // geometry
            mag_thickness_nm = 36.1;
            nonmag_thickness_nm = 21.32;
            theta = 0.; // external static magnetic field is assumed sufficiently strong that the magnetization aligns with the field

            h_static = 0.;
            // ac microwave intensity
            hpara2 = 2.;
            hperp2 = 2.;

            // sweep range
            min_B_Tesla = 1.3;
            max_B_Tesla = 2.6;
            n_B = 100;
            min_freq_GHz = 10.;
            max_freq_GHz = 50.;
            n_freq = 100;

            if (theta != 0.0) {
                throw std::runtime_error("Nonzero magnetization angle is not yet implemented");
            }
        }

        void setup_grids() {
            B_Tesla.linspace(min_B_Tesla, max_B_Tesla, n_B);
            freq_GHz.linspace(min_freq_GHz, max_freq_GHz, n_freq);
        }
    };

}