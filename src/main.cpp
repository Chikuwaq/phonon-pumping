#include <iostream>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <cmath>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#include "../include/constants.h"
#include "../include/vector.h"
#include "../include/matrix.h"
#include "../include/materialparameters.h"
#include "../include/geometry.h"

using namespace phonon_pumping;

static const std::string spectrum_datafile_name = "spectrum.dat";
static const std::string gnuplot_file_name = "plot_spectrum.plt";

void print_usage() {
    std::cout << "Usage: ./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file" << std::endl;
}

void calculate(const MaterialParameters& material, const Geometry& geometry, Matrix2D<double>& FMR_power_absorption) {
    // frequency arrays
    Vector<double> omega_theta_Hz; omega_theta_Hz.zeros(geometry.n_B);
    for (size_t i_B = 0; i_B < geometry.n_B; i_B++) {
        const double omH = material.omegaH(geometry.B_Tesla(i_B));
        const double omKMinusM = material.omegaK() - material.omegaM();
        omega_theta_Hz(i_B) = sqrt(pow(omH, 2) + pow(omKMinusM, 2) * cos(2 * geometry.theta) + (3. * pow(cos(geometry.theta), 2) - 1.) * omH * omKMinusM);
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

    Matrix2D<double> h; h.zeros(4, geometry.n_freq); // 4-element vector with elements being a function of frequency
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        h(0, i_freq) = -omega_Hz(i_freq) * (beta(1, 0, i_freq) * cos(phi(i_freq)) * cos(xi(i_freq)) - beta(0, 1, i_freq) * sin(phi(i_freq)) * sin(xi(i_freq))) - 0.5 * material.mag_eta_Hz * (beta_prime(0, 0, i_freq) * sin(phi(i_freq)) * cos(xi(i_freq)) + beta_prime(1, 1, i_freq) * cos(phi(i_freq)) * sin(xi(i_freq)));
        h(1, i_freq) = omega_Hz(i_freq) * (beta(0, 0, i_freq) * sin(phi(i_freq)) * cos(xi(i_freq)) + beta(1, 1, i_freq) * cos(phi(i_freq)) * sin(xi(i_freq))) - 0.5 * material.mag_eta_Hz * (beta_prime(1, 0, i_freq) * cos(phi(i_freq)) * cos(xi(i_freq)) - beta_prime(0, 1, i_freq) * sin(phi(i_freq)) * sin(xi(i_freq)));
        h(2, i_freq) = omega_Hz(i_freq) * ((cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta(0, 0, i_freq)) * cos(phi(i_freq)) * cos(xi(i_freq)) - (sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta(1, 1, i_freq)) * sin(phi(i_freq)) * sin(xi(i_freq)) - 2. * cosh(nonmag_damped_phase) * cos(xi(i_freq))) + 0.5 * material.mag_eta_Hz * ((sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta_prime(1, 0, i_freq)) * sin(phi(i_freq)) * cos(xi(i_freq)) + (cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta_prime(0, 1, i_freq)) * cos(phi(i_freq)) * sin(xi(i_freq)) - 2. * sinh(nonmag_damped_phase) * sin(xi(i_freq)));
        h(3, i_freq) = -omega_Hz(i_freq) * ((sinh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta(1, 0, i_freq)) * sin(phi(i_freq)) * cos(xi(i_freq)) + (cosh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta(0, 1, i_freq)) * cos(phi(i_freq)) * sin(xi(i_freq)) - 2. * sinh(nonmag_damped_phase) * sin(xi(i_freq))) + 0.5 * material.mag_eta_Hz * ((cosh(mag_damped_phase) * cosh(nonmag_damped_phase) + beta_prime(0, 0, i_freq)) * cos(phi(i_freq)) * cos(xi(i_freq)) - (sinh(mag_damped_phase) * sinh(nonmag_damped_phase) + beta_prime(1, 1, i_freq)) * sin(phi(i_freq)) * sin(xi(i_freq)) - 2. * cosh(nonmag_damped_phase) * cos(xi(i_freq)));
    }

    // dynamical magnetoelastic coupling 
    Vector<double> g_real; g_real.zeros(geometry.n_freq);
    Vector<double> g_imag; g_imag.zeros(geometry.n_freq);
    for (size_t i_freq = 0; i_freq < geometry.n_freq; i_freq++) {
        const double omega = omega_Hz(i_freq);
        const double g_front = -pow(material.omegac(), 2) * material.Ms_Jm3 / constants::GAMMA / mag_thickness_m / material.mag_rho_kgm3 / material.mag_ct_ms / (pow(omega, 2) + pow(0.5 * material.mag_eta_Hz, 2));
        g_real(i_freq) = g_front * (-omega * (h(0, i_freq) * h(3, i_freq) - h(1, i_freq) * h(2, i_freq)) + 0.5 * material.mag_eta_Hz * (h(0, i_freq) * h(2, i_freq) + h(1, i_freq) * h(3, i_freq)))/(pow(h(0, i_freq), 2) + pow(h(1, i_freq), 2));
        g_imag(i_freq) = g_front * (omega * (h(0, i_freq) * h(2, i_freq) + h(1, i_freq) * h(3, i_freq)) + 0.5 * material.mag_eta_Hz * (h(0, i_freq) * h(3, i_freq) - h(1, i_freq) * h(2, i_freq)))/(pow(h(0, i_freq), 2) + pow(h(1, i_freq), 2));
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
            Omega_real(0, i_freq, i_B) = omH + (omK-omM-g_real(i_freq)) * pow(cos(geometry.theta), 2);
            Omega_imag(0, i_freq, i_B) = -material.Gilbert_damping * omega_Hz(i_freq) - g_imag(i_freq) * pow(cos(geometry.theta), 2);
            Omega_real(1, i_freq, i_B) = omH + (omK-omM) * cos(2. * geometry.theta) - g_real(i_freq) * pow(cos(2 * geometry.theta), 2) - gl_real(i_freq) * pow(sin(2 * geometry.theta), 2);
            Omega_imag(1, i_freq, i_B) = -material.Gilbert_damping * omega_Hz(i_freq) - g_imag(i_freq) * pow(cos(2 * geometry.theta), 2) - gl_imag(i_freq) * pow(sin(2 * geometry.theta), 2);

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

void output_spectrum(const Vector<double>& x_values, const Vector<double>& y_values, const Matrix2D<double>& z_values) {
    std::cout << "\nWriting spectrum data to file..." << std::endl;

    std::ofstream file(spectrum_datafile_name);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open output file!");
    }
    assert(x_values.size() == z_values.nRow() && "Number of data points along x-axis does not match");
    assert(y_values.size() == z_values.nCol() && "Number of data points along y-axis does not match");
    
    for (size_t i_x = 0; i_x < x_values.size(); i_x++) {
        const double x = x_values(i_x);
        for (size_t i_y = 0; i_y < y_values.size(); i_y++) {
            const double y = y_values(i_y);
            file << std::setw(7) << x;
            file << std::setw(7) << y;
            file << std::setw(14) << z_values(i_x, i_y);
            file << std::endl;
        }
        file << std::endl; // gnuplot requires a blank line after every change in the x coordinate

    }
    file.close();
}

void write_gnuplot_script(const Geometry& geometry) {
    std::cout << "\nWriting gnuplot script file..." << std::endl;

    std::ofstream file(gnuplot_file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open output file!");
    }

    #if defined(__APPLE__)
    file << "set term pdfcairo enhanced color font 'Helvetica,18' size 5in,3in" << std::endl;
    #else
    file << "set term pdfcairo enhanced color font 'sans,18' size 5in,3in" << std::endl;
    #endif
    file << "set output 'FMR_MagThickness" << geometry.mag_thickness_nm << "_NonmagThickness" << geometry.nonmag_thickness_nm << ".pdf' " << std::endl;
    file << "set lmargin screen 0.20" << std::endl;
    file << "set rmargin screen 0.80" << std::endl;

    file << "set pm3d map" << std::endl;
    file << "set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000') # define color gradation. Relative to cbrange" << std::endl;

    file << "set title 'FMR power absorption: Magnet/Nonmagnet " << geometry.mag_thickness_nm << "/" << geometry.nonmag_thickness_nm << " nm'" << std::endl;

    file << "set xlabel 'Out-of-plane magnetic field [T]'" << std::endl;
    file << "set ylabel 'Frequency [GHz]'" << std::endl;

    file << "sp '" << spectrum_datafile_name << "' notitle w pm3d" << std::endl;

    for (auto& [label, value] : geometry.horizontal_lines) {
        file << "p " << value << " title '" << label << "' w line" << std::endl; // TODO: distinguish by colors
    }

    file << "unset pm3d" << std::endl;
    file << "unset term" << std::endl;
    file << "set output" << std::endl;
}

int main(int argc, char *argv[]) {

    if (argc < 5) {
        print_usage();
        return -1;
    }

    std::filesystem::path path_material_parameter;
    std::filesystem::path path_geometry;

    for(int count = 1; count < argc; count++ ) {
        std::string s = argv[count];
        if (s == "-m") {
            count++;
            path_material_parameter = argv[count];
        }
        else if (s == "-g") {
            count++;
            path_geometry = argv[count];
        }
    }

    if (path_material_parameter.empty()) {
        std::cout << "Path to material parameter file is empty!" << std::endl;
        print_usage();
        return -1;
    }
    if (path_geometry.empty()) {
        std::cout << "Path to geometry file is empty!" << std::endl;
        print_usage();
        return -1;
    }

    try {
        MaterialParameters material;
        material.set(path_material_parameter);

        Geometry geometry;
        geometry.load(path_geometry, material);

        // --- compute FMR absorption spectrum
        geometry.setup_grids();
        
        Matrix2D<double> FMR_power_absorption; 
        calculate(material, geometry, FMR_power_absorption);

        output_spectrum(geometry.B_Tesla, geometry.freq_GHz, FMR_power_absorption);
        write_gnuplot_script(geometry);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "\nPhononPumpingPlotter DONE. Run the gnuplot script " << gnuplot_file_name << " to generate figures." << std::endl;
    
    return 0;
}