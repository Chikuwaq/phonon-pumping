/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <iostream>
#include <filesystem>
#include <string>
#include <vector>

#include "parser.h"
#include "task.h"

using TaskID = phonon_pumping::task::TaskID;

std::vector<TaskID> phonon_pumping::task::load(const std::filesystem::path& filepath) {
    std::cout << "\nReading task file..." << std::endl;

    auto parsed_input = parse_file(filepath);
    std::vector<TaskID> result;

    if (parsed_input.find("plot_FMR_spectrum") != parsed_input.end() && parsed_input["plot_FMR_spectrum"] > 0) {
        result.push_back(TaskID::PLOT_FMR_SPECTRUM);
    }

    if (parsed_input.find("target_TA1_GHz") != parsed_input.end() && parsed_input["target_TA1_GHz"] > 0) {
        result.push_back(TaskID::FIT_TA_SOUND_VELOCITIES);
    }

    if (parsed_input.find("target_LA1_GHz") != parsed_input.end() && parsed_input["target_LA1_GHz"] > 0) {
        result.push_back(TaskID::FIT_LA_SOUND_VELOCITIES);
    }

    return result;
}

std::vector<double> phonon_pumping::task::extract_freqs(const std::filesystem::path& filepath, const enhance_mechanism mode) {
    std::string mode_label;
    if (mode == enhance_mechanism::EIGENMODE_TA) mode_label = "TA";
    else if (mode ==enhance_mechanism::EIGENMODE_LA) mode_label = "LA";

    constexpr u_int16_t max_num_frequencies = 10;

    std::cout << "\nReading task file..." << std::endl;

    auto parsed_input = parse_file(filepath);
    std::vector<double> frequencies;

    std::vector<std::string> keys;
    for (u_int16_t i = 1; i < max_num_frequencies; i++) {
        keys.push_back("target_" + mode_label + std::to_string(i) + "_GHz");
    }

    for (std::string& key : keys) {
        if (parsed_input.find(key) != parsed_input.end()) {
            frequencies.push_back(parsed_input[key]);
        }
    }

    return frequencies;
}