#include "gnuplot.h"


void phonon_pumping::output_spectrum(const Vector<double>& x_values, 
                                     const Vector<double>& y_values, 
                                     const Matrix2D<double>& z_values) {
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

double phonon_pumping::screen_position_x(const Geometry& geometry, 
                                         const double position_with_respect_to_screen) {
    return (1.0 - position_with_respect_to_screen) * geometry.min_B_Tesla + position_with_respect_to_screen * geometry.max_B_Tesla;
}