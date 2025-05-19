/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "materialparameters.h"

void phonon_pumping::MaterialParameters::set(const std::filesystem::path& filepath) {
	std::cout << "\nReading material parameter file..." << std::endl;

	auto result = parse_file(filepath);

	for (const auto& [key, value] : result) {
		std::cout << key << " = " << value << std::endl;
	}

	mag_rho_kgm3 = result["mag_rho_kgm3"];
	mag_ct_ms = result["mag_ct_ms"];
	mag_cl_ms = result["mag_cl_ms"];
	mag_eta_Hz = 2.0 * constants::PI * result["mag_eta_Hz"];

	nonmag_rho_kgm3 = result["nonmag_rho_kgm3"];
	nonmag_ct_ms = result["nonmag_ct_ms"];
	nonmag_cl_ms = result["nonmag_cl_ms"];
	nonmag_eta_Hz = 2.0 * constants::PI * result["nonmag_eta_Hz"];

	K1_Jm3 = result["K1_Jm3"];
	Ms_Jm3 = result["Ms_Jm3"];
	Gilbert_damping = result["Gilbert_damping"];

	b1_Jm3 = result["b1_Jm3"];
	b2_Jm3 = result["b2_Jm3"];
};