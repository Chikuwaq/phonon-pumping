/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "elasticity.h"
#include "numerics.h"

void phonon_pumping::calculate_Gilbert_damping_enhancing_freq(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw) {
    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann]) {
        std::cout << "\nSearching for frequencies at which the TA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        bool first_appearance = true;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            std::string label = "n_t=" + std::to_string(i);
            if (first_appearance) {
                label += " (N-N)";
                first_appearance = false;
            }
            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann].emplace_back(label, resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann]) {
        std::cout << "Searching for frequencies at which the LA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
        bool first_appearance = true;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_LA_both_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            std::string label = "n_l=" + std::to_string(i);
            if (first_appearance) {
                label += " (N-N)";
                first_appearance = false;
            }
            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann].emplace_back(label, resonance_freq_GHz);
        }
    }

    if (geometry.should_draw_frequencies[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann]) {
        std::cout << "Searching for frequencies at which the TA mode profile (Dirichlet and Neumann boundary conditions on each side) is in phase with the interface magnetoelastic stress ..." << std::endl;
        bool first_appearance = true;
        for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
            const double resonance_freq_GHz = material.stress_matching_freq_TA_Dirichlet_Neumann_GHz(geometry.mag_thickness_nm / constants::SCALE_TO_NANO, i);
            if (resonance_freq_GHz < geometry.min_freq_GHz) continue;
            if (geometry.max_freq_GHz < resonance_freq_GHz) break;

            std::string label = "n_t=" + std::to_string(i);
            if (first_appearance) {
                label += " (D-N)";
                first_appearance = false;
            }
            frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann].emplace_back(label, resonance_freq_GHz);
        }
    }
}

// analytical expression of the left-hand-side of Eq. (38) at given frequency
// Calling this method for FMR spectrum calculation is not recommended; calculation will take longer.
double phonon_pumping::eigenvalue_equation(const MaterialParameters& material, 
                                           const Geometry& geometry, 
                                           const bool is_TA,
                                           const double freq_GHz) {
    const double mag_sound_velocity = is_TA ? material.mag_ct_ms : material.mag_cl_ms;
    const double nonmag_sound_velocity = is_TA ? material.nonmag_ct_ms : material.nonmag_cl_ms;

    const double impedance_ratio = material.nonmag_rho_kgm3 * nonmag_sound_velocity / (material.mag_rho_kgm3 * mag_sound_velocity);
    const double omega_Hz = 2. * constants::PI * freq_GHz / constants::SCALE_TO_GIGA;
    const double phi = omega_Hz * (geometry.mag_thickness_nm / constants::SCALE_TO_NANO) / mag_sound_velocity; // dimensionless
    const double xi = omega_Hz * (geometry.nonmag_thickness_nm / constants::SCALE_TO_NANO) / nonmag_sound_velocity; // dimensionless
    return sin(phi) * cos(xi) + impedance_ratio * cos(phi) * sin(xi);
}

void phonon_pumping::find_acoustic_eigenfrequencies(const MaterialParameters& material, const Geometry& geometry, list_of_frequencies& frequencies_to_draw) {
    if (geometry.should_draw_frequencies[enhance_mechanism::EIGENMODE_TA]) {
        std::cout << "\nSearching for TA mode eigenfrequencies..." << std::endl;
        append_eigenfrequencies(material, geometry, enhance_mechanism::EIGENMODE_TA, frequencies_to_draw);
    }
    if (geometry.should_draw_frequencies[enhance_mechanism::EIGENMODE_LA]) {
        std::cout << "\nSearching for LA mode eigenfrequencies..." << std::endl;
        append_eigenfrequencies(material, geometry, enhance_mechanism::EIGENMODE_LA, frequencies_to_draw);
    }
}

void phonon_pumping::append_eigenfrequencies(const MaterialParameters& material, 
                                             const Geometry& geometry, 
                                             const enhance_mechanism i_mechanism, 
                                             list_of_frequencies& frequencies_to_draw) {
    // identify appropriate output label and eigenvalue equation
    std::string label;
    std::function<double(double)> lhs;
    if (i_mechanism == enhance_mechanism::EIGENMODE_TA) {
        label = "TA";
        lhs = [material, geometry](double freq_GHz) { return eigenvalue_equation(material, geometry, true, freq_GHz); };
    }
    else if (i_mechanism == enhance_mechanism::EIGENMODE_LA) {
        label = "LA";
        lhs = [material, geometry](double freq_GHz) { return eigenvalue_equation(material, geometry, false, freq_GHz); };
    }
    else {
        throw std::runtime_error("Gilbert-damping enhancement mechanism " + std::to_string(i_mechanism) + " is not supported");
    }
    
    // find roots of the eigenvalue equation by the Secant method and store them 
    const double dfreq = geometry.freq_GHz(1) - geometry.freq_GHz(0);
    uint16_t i_eigenmode = 0;
    for (double freq_GHz = 0; freq_GHz <= geometry.freq_GHz(geometry.n_freq - 1); freq_GHz += dfreq) { // loop must start from 0 GHz to enumerate eigenfrequencies
        const double freq0 = freq_GHz;
        const double freq1 = freq_GHz + dfreq;
        if (lhs(freq0) * lhs(freq1) > 0) continue;

        i_eigenmode++;
        if (freq_GHz < geometry.freq_GHz(0)) continue; // do not draw eigenfrequencies out of plot range

        const double eigenfrequency_GHz = numerics::secant_method(lhs, freq0, freq1);
        frequencies_to_draw[i_mechanism].emplace_back(label + " " + std::to_string(i_eigenmode), eigenfrequency_GHz);
    }
}