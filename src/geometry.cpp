#include <cmath>
#include <string>

#include "geometry.h"

uint16_t phonon_pumping::Geometry::expect_uint(const double value) {
	if (value < 0.0) {
		throw std::runtime_error("Negative value is not allowed for index!");
	}
	return static_cast<uint16_t>(std::floor(value));
}

// Range will be null if the key is absent in the geometry file.
std::tuple<uint16_t, uint16_t> phonon_pumping::Geometry::get_index_range(std::map<std::string, double>& parsed_input, const std::string& key) {
	uint16_t i_min = (parsed_input.find(key + "_min_index") != parsed_input.end()) ? expect_uint(parsed_input[key + "_min_index"]) : std::numeric_limits<uint16_t>::max();
	uint16_t i_max = (parsed_input.find(key + "_max_index") != parsed_input.end()) ? expect_uint(parsed_input[key + "_max_index"]) : std::numeric_limits<uint16_t>::max() - 1;
	return std::tuple(i_min, i_max);
}

void phonon_pumping::Geometry::load(const std::filesystem::path& filepath, const MaterialParameters& material) {
	std::cout << "\nReading geometry file..." << std::endl;

	auto parsed_input = parse_file(filepath);

	if (parsed_input.find("theta") == parsed_input.end()) {
		throw std::runtime_error("Magnetization angle is missing in the input file");
	}
	if (parsed_input["theta"] != 0.0) {
		throw std::runtime_error("Nonzero magnetization angle is not yet implemented");
	}
	for (const auto& [key, value] : parsed_input) {
		std::cout << key << " = " << value << std::endl;
	}

	mag_thickness_nm = parsed_input["mag_thickness_nm"];
	nonmag_thickness_nm = parsed_input["nonmag_thickness_nm"];
	theta = parsed_input["theta"];

	h_static = parsed_input["h_static"];
	hpara2 = parsed_input["hpara2"];
	hperp2 = parsed_input["hperp2"];

	min_B_Tesla = parsed_input["min_B_Tesla"];
	max_B_Tesla = parsed_input["max_B_Tesla"];
	n_B = parsed_input["n_B"];
	min_freq_GHz = parsed_input["min_freq_GHz"];
	max_freq_GHz = parsed_input["max_freq_GHz"];
	n_freq = parsed_input["n_freq"];

	// calculate resonant frequencies and their labels for horizontal line outputs
	std::cout << "Calculating resonant frequencies..." << std::endl;
	auto [i_min_TA_NN, i_max_TA_NN] = get_index_range(parsed_input, "stress_matching_freq_TA_both_Neumann");
	for (uint16_t i = i_min_TA_NN; i <= i_max_TA_NN; i++) {
		horizontal_lines.push_back(std::make_pair("Stress matching TA (Neumann) " + std::to_string(i), material.stress_matching_freq_TA_both_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i)));
	}

	auto [i_min_LA_NN, i_max_LA_NN] = get_index_range(parsed_input, "stress_matching_freq_LA_both_Neumann");
	for (uint16_t i = i_min_LA_NN; i <= i_max_LA_NN; i++) {
		horizontal_lines.push_back(std::make_pair("Stress matching LA (Neumann) " + std::to_string(i), material.stress_matching_freq_LA_both_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i)));
	}

	auto [i_min_TA_DN, i_max_TA_DN] = get_index_range(parsed_input, "stress_matching_freq_TA_Dirichlet_Neumann");
	for (uint16_t i = i_min_TA_DN; i <= i_max_TA_DN; i++) {
		horizontal_lines.push_back(std::make_pair("Stress matching TA (Dirichlet-Neumann) " + std::to_string(i), material.stress_matching_freq_TA_Dirichlet_Neumann_GHz(mag_thickness_nm / constants::SCALE_TO_NANO, i)));
	}
}

void phonon_pumping::Geometry::setup_grids() {
	B_Tesla.linspace(min_B_Tesla, max_B_Tesla, n_B);
	freq_GHz.linspace(min_freq_GHz, max_freq_GHz, n_freq);
}