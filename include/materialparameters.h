#include <filesystem>
#include <iostream>

#include "constants.h"
#include "parser.h"

namespace phonon_pumping {

    struct MaterialParameters {
        public:
        // elastic properties
        double mag_rho_kgm3;
        double mag_ct_ms;
        double mag_cl_ms;
        double mag_eta_Hz; // in general depends on freq_Hz, but approximately constant in the GHz range

        double nonmag_rho_kgm3;
        double nonmag_ct_ms;
        double nonmag_cl_ms;
        double nonmag_eta_Hz; // in general depends on freq_Hz, but approximately constant in the GHz range

        // magnetic properties
        double K1_Jm3;
        double Ms_Jm3;
        double Gilbert_damping;

        // magnetoelastic properties
        double b1_Jm3;
        double b2_Jm3;

        void set(const std::filesystem::path& filepath);

        inline double mag_kappat() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / mag_ct_ms;
        }

        inline double mag_kappal() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / mag_cl_ms;
        }

        inline double nonmag_kappat() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / nonmag_ct_ms;
        }

        inline double nonmag_kappal() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / nonmag_cl_ms;
        }

        inline double omegaK() const {
            return constants::GAMMA * 2.0 * K1_Jm3 / Ms_Jm3;
        }

        inline double omegaM() const {
            return constants::GAMMA * constants::MU0 * Ms_Jm3;
        }

        inline double omegac() const {
            return 0.5 * omegaM() + constants::GAMMA * (b2_Jm3 - K1_Jm3) / Ms_Jm3;
        }

        inline double omegacl() const {
            return constants::GAMMA * b1_Jm3 / Ms_Jm3;
        }

        inline double omegaH(const double magnetic_field) const {
            return constants::GAMMA * magnetic_field;
        }

        inline double stress_matching_freqt_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_ct_ms * (2. * mode - 1.) / 2. / thickness_m;
        }

        inline double stress_matching_freql_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_cl_ms * (2. * mode - 1.) / 2. / thickness_m;
        }

        inline double stress_matching_half_fixed_freqt_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_ct_ms * (2. * mode - 1.) / 4. / thickness_m;
        }
    };

}