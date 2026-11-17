/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "gnuplot.h"

void phonon_pumping::write_gnuplot_script(const MaterialParameters& material, 
                                          const Geometry& geometry, 
                                          const list_of_frequencies& frequencies_to_draw) {
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

    file << std::endl;
    file << "set title 'FMR power absorption: Magnet/Nonmagnet " << geometry.mag_thickness_nm << "/" << geometry.nonmag_thickness_nm << " nm'" << std::endl;

    file << "set xlabel 'Out-of-plane magnetic field [T]'" << std::endl;
    file << "set ylabel 'Frequency [GHz]'" << std::endl;
    
    file << "set xrange [" << geometry.min_B_Tesla << ":" << geometry.max_B_Tesla << "]" << std::endl;
    file << "set yrange [" << geometry.min_freq_GHz << ":" << geometry.max_freq_GHz << "]" << std::endl;

    write_gnuplot_command_draw_frequency_line(material, geometry, "yellow", frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann], false, file);
    write_gnuplot_command_draw_frequency_line(material, geometry, "olive", frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann], false, file);
    write_gnuplot_command_draw_frequency_line(material, geometry, "red", frequencies_to_draw[enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann], false, file);
    write_gnuplot_command_draw_frequency_line(material, geometry, "white", frequencies_to_draw[enhance_mechanism::EIGENMODE_TA], true, file);
    write_gnuplot_command_draw_frequency_line(material, geometry, "white", frequencies_to_draw[enhance_mechanism::EIGENMODE_LA], true, file);

    file << "sp '" << spectrum_datafile_name << "' notitle w pm3d" << std::endl;

    file << std::endl;
    file << "unset pm3d" << std::endl;
    file << "unset term" << std::endl;
    file << "set output" << std::endl;
}

void phonon_pumping::write_gnuplot_command_draw_frequency_line(const MaterialParameters& material, 
                                                               const Geometry& geometry, 
                                                               const std::string& color_name, 
                                                               const std::vector<std::pair<std::string, double>>& frequencies_to_draw, 
                                                               const bool outside, 
                                                               std::ofstream& file) {
    if (frequencies_to_draw.size() == 0) return;

    // estimate FMR frequency at the middle magnetic field
    const double middle_B_Tesla = 0.5 * (geometry.min_B_Tesla + geometry.max_B_Tesla);
    const double FMR_freq_at_middle_B_GHz = 0.5 * (material.omega11_no_MEC_no_Gilbert_damping(middle_B_Tesla, geometry.magnetization_angle_radian) + material.omega22_no_MEC_no_Gilbert_damping(middle_B_Tesla, geometry.magnetization_angle_radian)) * constants::SCALE_TO_GIGA / 2.0 / constants::PI;
    const double label_position_shifty = 0.05 * (geometry.max_freq_GHz - geometry.min_freq_GHz);

    constexpr double middle_space = 0.2;
    const double x0 = screen_position_x(geometry, 0.0);
    const double x1 = screen_position_x(geometry, (0.0 + 0.5 - middle_space / 2.0) / 2.0);
    const double x2 = screen_position_x(geometry, 0.5 - middle_space / 2.0);
    const double x3 = screen_position_x(geometry, 0.5 + middle_space / 2.0);
    const double x4 = screen_position_x(geometry, (0.5 + middle_space / 2.0 + 1.0) / 2.0);
    const double x5 = screen_position_x(geometry, 1.0);
    for (auto& [label, freq_GHz] : frequencies_to_draw) {
        // prevent the label from overlapping with the FMR peaks
        const bool use_right_half = (freq_GHz < FMR_freq_at_middle_B_GHz);
         
        double arrow_start, arrow_end;
        if (use_right_half) {
            arrow_start = outside ? x4 : x3;
            arrow_end = outside ? x5 : x4;
        }
        else {
            arrow_start = outside ? x0 : x1;
            arrow_end = outside ? x1 : x2;
        }
        file << "set arrow from " << arrow_start << "," << freq_GHz << " to " << arrow_end << "," << freq_GHz << " nohead front lc '" << color_name << "'" << std::endl;
        file << "set label '" << label << "' " << " at " << arrow_start << ", " << freq_GHz + label_position_shifty << " front font ',12' textcolor rgb '" << color_name << "'" << std::endl;
    }
}

double phonon_pumping::screen_position_x(
    const Geometry& geometry, 
    const double position_with_respect_to_screen
) {
    return (1.0 - position_with_respect_to_screen) * geometry.min_B_Tesla + position_with_respect_to_screen * geometry.max_B_Tesla;
}

// The data file uses gnuplot's "matrix nonuniform" layout:
//           x_0       x_1      ... x_{nRow-1}
//   y_0     M[0][0]   M[1][0]  ... M[nRow-1][0]
//   y_1     M[0][1]   M[1][1]  ... M[nRow-1][1]
void phonon_pumping::write_heatmap_xy(
    const std::vector<std::vector<double>>& data,
    const std::vector<double>& x, 
    const std::vector<double>& y,
    const std::string& x_label,
    const std::string& y_label,
    const std::string& cb_label,
    const std::string& data_file,
    const std::string& script_file,
    const std::string& image_file
) {
    const size_t nCol = data[0].size(); // = number of y values
    const size_t nRow = data.size(); // = number of x values

    // ---- data file ----
    std::ofstream dat(data_file);
    if (!dat) throw std::runtime_error("Cannot open file for writing: " + data_file);

    dat << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);

    dat << "Result";
    for (const double x_value : x) {
        dat << " " << x_value;
    }
    dat << "\n";

    for (size_t i_y = 0; i_y < nCol; i_y++) {
        dat << y[i_y];
        for (size_t i_x = 0; i_x < nRow; i_x++) {
            dat << " " << data[i_x][i_y];
        }
        dat << "\n";
    }
    dat.close();
    if (!dat) throw std::runtime_error("Failed to write file: " + data_file);

    // ---- gnuplot script ----
    std::ofstream gp(script_file);
    if (!gp) throw std::runtime_error("Cannot open file for writing: " + script_file);

    gp << "set terminal pngcairo size 800,700\n"
        << "set output '" << image_file << "'\n"
        << "set xlabel '" << x_label << "'\n"
        << "set ylabel '" << y_label << "'\n"
        << "set cblabel '" << cb_label << "'\n"
        << gnuplot_palette_cividis()
        << "set autoscale fix\n"
        << "plot '" << data_file << "' nonuniform matrix using 1:2:3 with image notitle\n";
    gp.close();
    if (!gp) throw std::runtime_error("Failed to write file: " + script_file);
}