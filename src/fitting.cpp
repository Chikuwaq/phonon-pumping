#include "fitting.h"
#include "console.h"


void phonon_pumping::fit_sound_velocities(
    const MaterialParameters& material, 
    const Geometry& geometry, 
    const AcousticFrequencyTable& target_freq_table, 
    const enhance_mechanism i_mechanism
) {
    const bool need_TA = (i_mechanism == enhance_mechanism::EIGENMODE_TA);
    const bool need_LA = (i_mechanism == enhance_mechanism::EIGENMODE_LA);
    const std::string mode_label = get_mode_label(i_mechanism);

    // copy data not to damage the database
    MaterialParameters materialTemp = material;
    double& sound_velocity_magnet = need_TA ? materialTemp.mag_ct_ms : materialTemp.mag_cl_ms;
    double& sound_velocity_nonmagnet = need_TA ? materialTemp.nonmag_ct_ms : materialTemp.nonmag_cl_ms;

    Geometry geometry_temp = geometry;
    
    const double min_sound_velocity = need_TA ? material.min_transverse_sound_velocity_ms : material.min_longitudinal_sound_velocity_ms;
    const double max_sound_velocity = need_TA ? material.max_transverse_sound_velocity_ms : material.max_longitudinal_sound_velocity_ms;

    constexpr double delta_sound_velocities = 10.0; // tentative
    const double n_sound_velocities_candidate = (max_sound_velocity - min_sound_velocity) / delta_sound_velocities;
    if (n_sound_velocities_candidate - floor(n_sound_velocities_candidate) > 0.1) {
        throw std::runtime_error("The range [min_sound_velocity, max_sound_velocity] must be a multiple of " + std::to_string(delta_sound_velocities));
    }
    const size_t n_sound_velocities = static_cast<size_t>(floor((max_sound_velocity - min_sound_velocity) / delta_sound_velocities));

    Vector<double> sound_velocities_magnet; 
    sound_velocities_magnet.linspace(min_sound_velocity, max_sound_velocity, n_sound_velocities);
    
    Vector<double> sound_velocities_nonmagnet = sound_velocities_magnet;

    Matrix2D<double> deviations;
    deviations.zeros(n_sound_velocities, n_sound_velocities);

    for (auto& [bilayer, target_freqs_for_one_bilayer_design] : target_freq_table) {
        geometry_temp.bilayer = bilayer;
        const std::vector<std::pair<std::string, double>>& target = target_freqs_for_one_bilayer_design[i_mechanism];
        const size_t n_targets = target.size();
        if (n_targets == 0) continue;

        for (uint16_t i_sound_velocity_magnet = 0; i_sound_velocity_magnet < n_sound_velocities; i_sound_velocity_magnet++) {
            sound_velocity_magnet = sound_velocities_magnet(i_sound_velocity_magnet);
    
            for (uint16_t i_sound_velocity_nonmagnet = 0; i_sound_velocity_nonmagnet < n_sound_velocities; i_sound_velocity_nonmagnet++) {
                sound_velocity_nonmagnet = sound_velocities_nonmagnet(i_sound_velocity_nonmagnet);
    
                list_of_frequencies eigenfrequencies;
                find_acoustic_eigenfrequencies(materialTemp, geometry_temp, need_TA, need_LA, eigenfrequencies);
    
                std::vector<double> mode_frequencies = extract_mode_frequencies(i_mechanism, eigenfrequencies);
                if (mode_frequencies.size() < n_targets) {
                    const double highest_freq_found = mode_frequencies[mode_frequencies.size() - 1];
                    const auto [mode_key, highest_freq_target] = target[n_targets - 1];
                    const bool should_find_next_freq = (highest_freq_found < highest_freq_target);
                    if (should_find_next_freq)
                        throw std::runtime_error("Number of eigenfrequencies found is smaller than the number of target frequencies! Please consider increasing 'max_freq_GHz' in the geometry file.");
                }
    
                const uint16_t n_freqs = std::min(n_targets, mode_frequencies.size());
                double sum = 0.0;
                for (uint16_t i_freq = 0; i_freq < n_freqs; i_freq++) {
                    const auto [mode_key, target_freq] = target[i_freq];
                    const double deviation = target_freq - mode_frequencies[i_freq];
                    sum += deviation * deviation;
                }
                deviations(i_sound_velocity_magnet, i_sound_velocity_nonmagnet) += std::sqrt(sum);
                // std::cout << "deviation(v_m = " << i_sound_velocity_magnet << ", v_nm = " << i_sound_velocity_nonmagnet << ") = " << std::sqrt(sum) << std::endl;
            }
        }
    }

    // report results
    auto [irow, icol, min_value] = deviations.min();
    std::cout << "Fitting results ";
    if (need_TA) std::cout << "TA modes:" << std::endl;
    else if (need_LA) std::cout << "LA modes:" << std::endl;

    std::cout << "\tSound velocity of magnet = " << sound_velocities_magnet(irow) << std::endl;
    std::cout << "\tSound velocity of nonmagnet = " << sound_velocities_nonmagnet(icol) << std::endl;
    std::cout << "\tsqrt(Sum of squared frequency deviation) = " << min_value << std::endl;

    // write deviations to file
    deviations.write_gnuplot_xy(
        sound_velocities_magnet,
        sound_velocities_nonmagnet,
        "Magnet sound velocity (m/s)", 
        "Nonmagnet sound velocity (m/s)", 
        "sqrt(sum of frequency deviations) (GHz)",
        "fitting_result_" + mode_label + "_sound_velocity.dat", 
        "plot_" + mode_label + "_fitting_deviations.plt", 
        mode_label + "_fitting_deviations.png"
    );

    wait_for_enter(); // the user might want to note down
}