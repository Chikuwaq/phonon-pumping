#include "elasticity.h"

void phonon_pumping::calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw) {
    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann]) {
        std::cout << "Searching frequencies at which the TA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann].emplace_back("n_t=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann]) {
        std::cout << "Searching frequencies at which the LA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_LA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann].emplace_back("n_l=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann]) {
        std::cout << "Searching frequencies at which the TA mode profile (Dirichlet and Neumann boundary conditions on each side) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_Dirichlet_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann].emplace_back("n_t=" + std::to_string(i) + " (Dirichlet-Neumann)", resonance_freq_GHz);
        }
    }
}