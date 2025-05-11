#include <filesystem>

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

        void set(const std::filesystem::path& filepath) {
            std::cout << "\nReading material parameter file..." << std::endl;

            auto result = parse_file(filepath);

            for (const auto& [key, value] : result) {
                std::cout << key << " = " << value << std::endl;
            }

            mag_rho_kgm3 = result["mag_rho_kgm3"];
            mag_ct_ms = result["mag_ct_ms"];
            mag_cl_ms = result["mag_cl_ms"];
            mag_eta_Hz = result["mag_eta_Hz"];

            nonmag_rho_kgm3 = result["nonmag_rho_kgm3"];
            nonmag_ct_ms = result["nonmag_ct_ms"];
            nonmag_cl_ms = result["nonmag_cl_ms"];
            nonmag_eta_Hz = result["nonmag_eta_Hz"];

            K1_Jm3 = result["K1_Jm3"];
            Ms_Jm3 = result["Ms_Jm3"];
            Gilbert_damping = result["Gilbert_damping"];

            b1_Jm3 = result["b1_Jm3"];
            b2_Jm3 = result["b2_Jm3"];
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