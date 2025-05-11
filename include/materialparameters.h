#include "constants.h"

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

        void set() {
            mag_rho_kgm3 = 8.8e3;
            mag_ct_ms = 1900.0; // to be verified
            mag_cl_ms = 2000.0; // not relevant for perpendicular magnetization
            mag_eta_Hz = 2.0 * constants::PI * 1e9;
        
            nonmag_rho_kgm3 = 2.145e4;
            nonmag_ct_ms = 1680.0; // to be verified
            nonmag_cl_ms = 2000.0; // not relevant for perpendicular magnetization
            nonmag_eta_Hz = 2.0 * constants::PI * 1e10;
        
            // magnetic properties
            K1_Jm3 = 0.513e6;
            Ms_Jm3 = 1400e3;
            Gilbert_damping = 5.8e-3;
        
            // magnetoelastic properties
            b1_Jm3 = -8.1e6;  // temporarily using B1 of hexiagonal crystal
            b2_Jm3 = 37.4e6;
        };

        double mag_kappat() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / mag_ct_ms;
        }

        double mag_kappal() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / mag_cl_ms;
        }

        double nonmag_kappat() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / nonmag_ct_ms;
        }

        double nonmag_kappal() const { // in general depends on freq_Hz, but approximately constant in the GHz range
            return mag_eta_Hz / 2.0 / nonmag_cl_ms;
        }

        double omegaK() const {
            return constants::GAMMA * 2.0 * K1_Jm3 / Ms_Jm3;
        }

        double omegaM() const {
            return constants::GAMMA * constants::MU0 * Ms_Jm3;
        }

        double omegac() const {
            return 0.5 * omegaM() + constants::GAMMA * (b2_Jm3 - K1_Jm3) / Ms_Jm3;
        }

        double omegacl() const {
            return constants::GAMMA * b1_Jm3 / Ms_Jm3;
        }

        double omegaH(const double magnetic_field) const {
            return constants::GAMMA * magnetic_field;
        }

        double stress_matching_freqt_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_ct_ms * (2. * mode - 1.) / 2. / thickness_m;
        }

        double stress_matching_freql_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_cl_ms * (2. * mode - 1.) / 2. / thickness_m;
        }

        double stress_matching_half_fixed_freqt_GHz(const double thickness_m, const uint16_t mode) {
            return constants::SCALE_TO_GIGA * mag_ct_ms * (2. * mode - 1.) / 4. / thickness_m;
        }
    };

}