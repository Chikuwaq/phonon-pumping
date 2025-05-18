
#include "magnetoelasticity.h"

void phonon_pumping::calculate_spectrum(const MaterialParameters& material, const Geometry& geometry, Matrix2D<double>& FMR_power_absorption) {
    // frequency arrays
    Vector<double> omega_theta_Hz; omega_theta_Hz.zeros(geometry.n_B);
    for (size_t i_B = 0; i_B < geometry.n_B; i_B++) {
        const double omH = material.omegaH(geometry.B_Tesla(i_B));
        const double omKMinusM = material.omegaK() - material.omegaM();
        omega_theta_Hz(i_B) = sqrt(pow(omH, 2) + pow(omKMinusM, 2) * cos(2 * geometry.magnetization_angle_radian) + (3. * pow(cos(geometry.magnetization_angle_radian), 2) - 1.) * omH * omKMinusM);
    }

    Vector<double> omega_Hz; omega_Hz.zeros(geometry.n_freq);
    Vector<double> phi; phi.zeros(geometry.n_freq);
    Vector<double> xi; xi.zeros(geometry.n_freq);
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        omega_Hz(i_freq) = 2. * constants::PI * geometry.freq_GHz(i_freq) / constants::SCALE_TO_GIGA;
        phi(i_freq) = omega_Hz(i_freq) * (geometry.mag_thickness_nm / constants::SCALE_TO_NANO) / material.mag_ct_ms; // dimensionless
        xi(i_freq) = omega_Hz(i_freq) * (geometry.nonmag_thickness_nm / constants::SCALE_TO_NANO) / material.nonmag_ct_ms; // dimensionless
    }

    const double impedance_ratio_t = material.nonmag_rho_kgm3 * material.nonmag_ct_ms / material.mag_rho_kgm3 / material.mag_ct_ms;
    const double mag_thickness_m = geometry.mag_thickness_nm / constants::SCALE_TO_NANO;
    const double nonmag_thickness_m = geometry.nonmag_thickness_nm / constants::SCALE_TO_NANO;
    const double mag_damped_phase = material.mag_kappat() * mag_thickness_m;
    const double nonmag_damped_phase = material.nonmag_kappat() * nonmag_thickness_m;

    Matrix3D<double> beta; beta.zeros(2, 2, geometry.n_freq); // 2x2 matrix with elements being a function of frequency
    Matrix3D<double> beta_prime; beta_prime.zeros(2, 2, geometry.n_freq); // 2x2 matrix with elements being a function of frequency
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        const double f_Hz = geometry.freq_GHz(i_freq) / constants::SCALE_TO_GIGA;
        beta(0, 0, i_freq) = cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + impedance_ratio_t * sinh(mag_damped_phase) * sinh(nonmag_damped_phase);
        beta(0, 1, i_freq) = cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + impedance_ratio_t * sinh(mag_damped_phase) * cosh(nonmag_damped_phase);
        beta(1, 0, i_freq) = sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + impedance_ratio_t * cosh(mag_damped_phase) * sinh(nonmag_damped_phase);
        beta(1, 1, i_freq) = sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + impedance_ratio_t * cosh(mag_damped_phase) * cosh(nonmag_damped_phase);
        beta_prime(0, 0, i_freq) = cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + impedance_ratio_t * material.nonmag_eta_Hz / material.mag_eta_Hz * sinh(mag_damped_phase) * sinh(nonmag_damped_phase);
        beta_prime(0, 1, i_freq) = cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + impedance_ratio_t * material.nonmag_eta_Hz / material.mag_eta_Hz * sinh(mag_damped_phase) * cosh(nonmag_damped_phase);
        beta_prime(1, 0, i_freq) = sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + impedance_ratio_t * material.nonmag_eta_Hz / material.mag_eta_Hz * cosh(mag_damped_phase) * sinh(nonmag_damped_phase);
        beta_prime(1, 1, i_freq) = sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + impedance_ratio_t * material.nonmag_eta_Hz / material.mag_eta_Hz * cosh(mag_damped_phase) * cosh(nonmag_damped_phase);
    }

    Matrix2D<double> ht; ht.zeros(4, geometry.n_freq); // 4-element vector with elements being a function of frequency
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        ht(0, i_freq) = -omega_Hz(i_freq) * (beta(1, 0, i_freq) * cos(phi(i_freq)) * cos(xi(i_freq)) - beta(0, 1, i_freq) * sin(phi(i_freq)) * sin(xi(i_freq))) - 0.5 * material.mag_eta_Hz * (beta_prime(0, 0, i_freq) * sin(phi(i_freq)) * cos(xi(i_freq)) + beta_prime(1, 1, i_freq) * cos(phi(i_freq)) * sin(xi(i_freq)));
        ht(1, i_freq) = omega_Hz(i_freq) * (beta(0, 0, i_freq) * sin(phi(i_freq)) * cos(xi(i_freq)) + beta(1, 1, i_freq) * cos(phi(i_freq)) * sin(xi(i_freq))) - 0.5 * material.mag_eta_Hz * (beta_prime(1, 0, i_freq) * cos(phi(i_freq)) * cos(xi(i_freq)) - beta_prime(0, 1, i_freq) * sin(phi(i_freq)) * sin(xi(i_freq)));
        ht(2, i_freq) = omega_Hz(i_freq) * ((cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta(0, 0, i_freq)) * cos(phi(i_freq)) * cos(xi(i_freq)) - (sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta(1, 1, i_freq)) * sin(phi(i_freq)) * sin(xi(i_freq)) - 2. * cosh(nonmag_damped_phase) * cos(xi(i_freq))) + 0.5 * material.mag_eta_Hz * ((sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta_prime(1, 0, i_freq)) * sin(phi(i_freq)) * cos(xi(i_freq)) + (cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta_prime(0, 1, i_freq)) * cos(phi(i_freq)) * sin(xi(i_freq)) - 2. * sinh(nonmag_damped_phase) * sin(xi(i_freq)));
        ht(3, i_freq) = -omega_Hz(i_freq) * ((sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta(1, 0, i_freq)) * sin(phi(i_freq)) * cos(xi(i_freq)) + (cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta(0, 1, i_freq)) * cos(phi(i_freq)) * sin(xi(i_freq)) - 2. * sinh(nonmag_damped_phase) * sin(xi(i_freq))) + 0.5 * material.mag_eta_Hz * ((cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta_prime(0, 0, i_freq)) * cos(phi(i_freq)) * cos(xi(i_freq)) - (sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta_prime(1, 1, i_freq)) * sin(phi(i_freq)) * sin(xi(i_freq)) - 2. * cosh(nonmag_damped_phase) * cos(xi(i_freq)));
    }

    // dynamical magnetoelastic coupling 
    Vector<double> gt_real; gt_real.zeros(geometry.n_freq);
    Vector<double> gt_imag; gt_imag.zeros(geometry.n_freq);
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        const double omega = omega_Hz(i_freq);
        const double gt_front = -pow(material.omegact(), 2) * material.Ms_Jm3 / constants::GAMMA / mag_thickness_m / material.mag_rho_kgm3 / material.mag_ct_ms / (pow(omega, 2) + pow(0.5 * material.mag_eta_Hz, 2));
        gt_real(i_freq) = gt_front * (-omega * (ht(0, i_freq) * ht(3, i_freq) - ht(1, i_freq) * ht(2, i_freq)) + 0.5 * material.mag_eta_Hz * (ht(0, i_freq) * ht(2, i_freq) + ht(1, i_freq) * ht(3, i_freq)))/(pow(ht(0, i_freq), 2) + pow(ht(1, i_freq), 2));
        gt_imag(i_freq) = gt_front * (omega * (ht(0, i_freq) * ht(2, i_freq) + ht(1, i_freq) * ht(3, i_freq)) + 0.5 * material.mag_eta_Hz * (ht(0, i_freq) * ht(3, i_freq) - ht(1, i_freq) * ht(2, i_freq)))/(pow(ht(0, i_freq), 2) + pow(ht(1, i_freq), 2));
    }
    Vector<double> gl_real; gl_real.zeros(geometry.n_freq);
    Vector<double> gl_imag; gl_imag.zeros(geometry.n_freq);
    // TODO: implement longitudinal mode and coupling
    // for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
    //     let omega = omega_Hz(i_freq);
    //     let g_front = -material.omegac().powf(2.) * material.Ms_Jm3 / constants::GAMMA / mag_thickness_m / material.mag_rho_kgm3 / material.mag_ct_ms / (pow(omega, 2) + (0.5 * material.mag_eta_Hz).powf(2.));
    //     gl_real(i_freq) = g_front * (-omega * (h(0, i_freq) * h(3, i_freq) - h(1, i_freq) * h(2, i_freq)) + 0.5 * material.mag_eta_Hz * (h(0, i_freq) * h(2, i_freq) + h(1, i_freq) * h(3, i_freq)))/(pow(h(0, i_freq), 2) + pow(h(1, i_freq), 2));
    //     gl_imag(i_freq) = g_front * (omega * (h(0, i_freq) * h(2, i_freq) + h(1, i_freq) * h(3, i_freq)) + 0.5 * material.mag_eta_Hz * (h(0, i_freq) * h(3, i_freq) - h(1, i_freq) * h(2, i_freq)))/(pow(h(0, i_freq), 2) + pow(h(1, i_freq), 2));
    // }

    // FMR spectrum
    Matrix3D<double> Omega_real; Omega_real.zeros(2, geometry.n_freq, geometry.n_B);
    Matrix3D<double> Omega_imag; Omega_imag.zeros(2, geometry.n_freq, geometry.n_B);
    Matrix2D<double> chi_tot_real; chi_tot_real.zeros(2, 2);
    Matrix2D<double> chi_tot_imag; chi_tot_imag.zeros(2, 2);
    FMR_power_absorption.zeros(geometry.n_B, geometry.n_freq);

    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        for (size_t i_B = 0; i_B < geometry.n_B; i_B++) {
            const double omH = material.omegaH(geometry.B_Tesla(i_B));
            const double omK = material.omegaK();
            const double omM = material.omegaM();
            Omega_real(0, i_freq, i_B) = omH + (omK-omM-gt_real(i_freq)) * pow(cos(geometry.magnetization_angle_radian), 2);
            Omega_imag(0, i_freq, i_B) = -material.Gilbert_damping * omega_Hz(i_freq) - gt_imag(i_freq) * pow(cos(geometry.magnetization_angle_radian), 2);
            Omega_real(1, i_freq, i_B) = omH + (omK-omM) * cos(2. * geometry.magnetization_angle_radian) - gt_real(i_freq) * pow(cos(2 * geometry.magnetization_angle_radian), 2) - gl_real(i_freq) * pow(sin(2 * geometry.magnetization_angle_radian), 2);
            Omega_imag(1, i_freq, i_B) = -material.Gilbert_damping * omega_Hz(i_freq) - gt_imag(i_freq) * pow(cos(2 * geometry.magnetization_angle_radian), 2) - gl_imag(i_freq) * pow(sin(2 * geometry.magnetization_angle_radian), 2);

            const double Delta_real = Omega_real(0, i_freq, i_B) * Omega_real(1, i_freq, i_B) - Omega_imag(0, i_freq, i_B) * Omega_imag(1, i_freq, i_B) - pow(omega_Hz(i_freq), 2);
            const double Delta_imag = Omega_real(0, i_freq, i_B) * Omega_imag(1, i_freq, i_B) + Omega_imag(0, i_freq, i_B) * Omega_real(1, i_freq, i_B);
            const double Delta_squared = pow(Delta_real, 2) + pow(Delta_imag, 2);
            chi_tot_real(0, 0) = constants::GAMMA * constants::MU0 * (Omega_real(0, i_freq, i_B) * Delta_real + Omega_imag(0, i_freq, i_B) * Delta_imag) / Delta_squared;
            chi_tot_imag(0, 0) = constants::GAMMA * constants::MU0 * (-Omega_real(0, i_freq, i_B) * Delta_imag + Omega_imag(0, i_freq, i_B) * Delta_real) / Delta_squared;
            chi_tot_real(0, 1) = constants::GAMMA * constants::MU0 * (-omega_Hz(i_freq)*Delta_imag) / Delta_squared;
            chi_tot_imag(0, 1) = constants::GAMMA * constants::MU0 * (-omega_Hz(i_freq)*Delta_real) / Delta_squared;
            chi_tot_real(1, 0) = -chi_tot_real(0, 1);
            chi_tot_imag(1, 0) = -chi_tot_imag(0, 1);
            chi_tot_real(1, 1) = constants::GAMMA * constants::MU0 * (Omega_real(1, i_freq, i_B) * Delta_real + Omega_imag(1, i_freq, i_B) * Delta_imag) / Delta_squared;
            chi_tot_imag(1, 1) = constants::GAMMA * constants::MU0 * (-Omega_real(1, i_freq, i_B) * Delta_imag + Omega_imag(1, i_freq, i_B) * Delta_real) / Delta_squared;

            FMR_power_absorption(i_B, i_freq) = chi_tot_imag(0, 0) * geometry.hpara2 * (geometry.hpara2 + geometry.h_static) + chi_tot_imag(1, 1) * geometry.hperp2 * geometry.hperp2 - chi_tot_imag(0, 1) * geometry.hperp2 * geometry.h_static;
        }
    }
}

