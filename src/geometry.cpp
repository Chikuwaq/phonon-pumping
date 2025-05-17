#include <cmath>
#include <string>

#include "geometry.h"

uint16_t phonon_pumping::Geometry::expect_uint(const double value) {
	if (value < 0.0) {
		throw std::runtime_error("Negative value is not allowed for index!");
	}
	return static_cast<uint16_t>(std::floor(value));
}

void phonon_pumping::Geometry::load(const std::filesystem::path& filepath, const MaterialParameters& material) {
	std::cout << "\nReading geometry file..." << std::endl;

	auto parsed_input = parse_file(filepath);

	if (parsed_input.find("magnetization_angle_radian") == parsed_input.end()) {
		throw std::runtime_error("Magnetization angle is missing in the input file");
	}
	if (parsed_input["magnetization_angle_radian"] != 0.0) {
		throw std::runtime_error("Nonzero magnetization angle is not yet implemented");
	}
	for (const auto& [key, value] : parsed_input) {
		std::cout << key << " = " << value << std::endl;
	}

	mag_thickness_nm = parsed_input["mag_thickness_nm"];
	nonmag_thickness_nm = parsed_input["nonmag_thickness_nm"];
	magnetization_angle_radian = parsed_input["magnetization_angle_radian"];

	h_static = parsed_input["h_static"];
	hpara2 = parsed_input["hpara2"];
	hperp2 = parsed_input["hperp2"];

	min_B_Tesla = parsed_input["min_B_Tesla"];
	max_B_Tesla = parsed_input["max_B_Tesla"];
	n_B = parsed_input["n_B"];
	min_freq_GHz = parsed_input["min_freq_GHz"];
	max_freq_GHz = parsed_input["max_freq_GHz"];
	n_freq = parsed_input["n_freq"];

	// calculate resonant frequencies in the output frequency range and store with their labels
	// for horizontal line outputs
	if (parsed_input.find("stress_matching_freq_TA_both_Neumann") != parsed_input.end()) {
		if (parsed_input["stress_matching_freq_TA_both_Neumann"] > 0) {
			std::cout << "Searching frequencies at which the TA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
			for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
				const double resonance_freq_GHz = material.stress_matching_freq_TA_both_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i);
				if (resonance_freq_GHz < min_freq_GHz) continue;
				if (max_freq_GHz < resonance_freq_GHz) break;

				horizontal_lines[Enhance_mechanism::STRESS_MATCHING_TA_NeumannNeumann].emplace_back("n_t=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
			}
		}
	}

	if (parsed_input.find("stress_matching_freq_LA_both_Neumann") != parsed_input.end()) {
		if (parsed_input["stress_matching_freq_LA_both_Neumann"] > 0) {
			std::cout << "Searching frequencies at which the LA mode profile (Neumann boundary condition on both sides) is in phase with the interface magnetoelastic stress ..." << std::endl;
			for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
				const double resonance_freq_GHz = material.stress_matching_freq_LA_both_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i);
				if (resonance_freq_GHz < min_freq_GHz) continue;
				if (max_freq_GHz < resonance_freq_GHz) break;

				horizontal_lines[Enhance_mechanism::STRESS_MATCHING_LA_NeumannNeumann].emplace_back("n_l=" + std::to_string(i) + " (Neumann-Neumann)", resonance_freq_GHz);
			}
		}
	}

	if (parsed_input.find("stress_matching_freq_TA_Dirichlet_Neumann") != parsed_input.end()) {
		if (parsed_input["stress_matching_freq_TA_Dirichlet_Neumann"] > 0) {
			std::cout << "Searching frequencies at which the TA mode profile (Dirichlet and Neumann boundary conditions on each side) is in phase with the interface magnetoelastic stress ..." << std::endl;
			for (uint16_t i = 1; i <= std::numeric_limits<uint16_t>::max(); i++) {
				const double resonance_freq_GHz = material.stress_matching_freq_TA_Dirichlet_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i);
				if (resonance_freq_GHz < min_freq_GHz) continue;
				if (max_freq_GHz < resonance_freq_GHz) break;

				horizontal_lines[Enhance_mechanism::STRESS_MATCHING_TA_DirichletNeumann].emplace_back("n_t=" + std::to_string(i) + " (Dirichlet-Neumann)", resonance_freq_GHz);
			}
		}
	}
}

void phonon_pumping::Geometry::setup_grids() {
	B_Tesla.linspace(min_B_Tesla, max_B_Tesla, n_B);
	freq_GHz.linspace(min_freq_GHz, max_freq_GHz, n_freq);
}