#include "elasticity.h"
#include "numerics.h"

void phonon_pumping::calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw) {
    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann]) {
        std::cout << "\nSearching for frequencies at which the TA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann].emplace_back("n_t=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann]) {
        std::cout << "Searching for frequencies at which the LA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_LA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann].emplace_back("n_l=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann]) {
        std::cout << "Searching for frequencies at which the TA mode profile (Dirichlet and Neumann boundary conditions on each side) is in phase with the interface magnetoelastic stress ..." << std::endl;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_Dirichlet_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann].emplace_back("n_t=" + std::to_string(i) + " (Dirichlet-Neumann)", resonance_freq_GHz);
        }
    }
}

// analytical expression of the left-hand-side of Eq. (38)
// Calling this method for FMR spectrum calculation is not recommended; calculation will take longer.
double phonon_pumping::lhs_Eq38(const MaterialParameters& material, const Geometry& geometry, const double freq_GHz) {
    const double impedance_ratio = material.nonmag_rho_kgm3 * material.nonmag_ct_ms / (material.mag_rho_kgm3 * material.mag_ct_ms);
    const double omega_Hz = 2. * constants::PI * freq_GHz / constants::SCALE_TO_GIGA;
    const double phi = omega_Hz * (geometry.mag_thickness_nm / constants::SCALE_TO_NANO) / material.mag_ct_ms; // dimensionless
    const double xi = omega_Hz * (geometry.nonmag_thickness_nm / constants::SCALE_TO_NANO) / material.nonmag_ct_ms; // dimensionless
    return sin(phi) * cos(xi) + impedance_ratio * cos(phi) * sin(xi);
}

void phonon_pumping::find_acoustic_eigenfrequencies(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw) {
    std::cout << "\nSearching for TA mode eigenfrequencies..." << std::endl;
    
    const double dfreq = geometry.freq_GHz(1) - geometry.freq_GHz(0);
    uint16_t i_eigenmode = 0;
    for (double freq_GHz = 0; freq_GHz <= geometry.freq_GHz(geometry.n_freq - 1); freq_GHz += dfreq) { // loop must start from 0 GHz to enumerate eigenfrequencies
        const double freq0 = freq_GHz;
        const double freq1 = freq_GHz + dfreq;
        if (lhs_Eq38(material, geometry, freq0) * lhs_Eq38(material, geometry, freq1) > 0) continue;

        i_eigenmode++;
        if (freq_GHz < geometry.freq_GHz(0)) continue; // do not draw eigenfrequencies out of plot range

        auto func = [material, geometry](double freq_GHz) { return lhs_Eq38(material, geometry, freq_GHz); };
        const double eigenfrequency_GHz = numerics::secant_method(func, freq0, freq1);
        frequencies_to_draw[enhance_mechanism::EIGENMODE_TA].emplace_back("TA mode " + std::to_string(i_eigenmode), eigenfrequency_GHz);
    }

    // TODO: find LA mode eigenfrequencies
}