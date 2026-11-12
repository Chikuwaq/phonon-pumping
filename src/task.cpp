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

    if (parsed_input.find("target_LA1_GHz") != parsed_input.end() && parsed_input["target_LA1_GHz"] > 0) {
        result.push_back(TaskID::FIT_SOUND_VELOCITIES);
    }

    return result;
}

std::vector<double> phonon_pumping::task::extract_freqs(const std::filesystem::path& filepath) {
    std::cout << "\nReading task file..." << std::endl;

    auto parsed_input = parse_file(filepath);
    std::vector<double> frequencies;

    std::vector<std::string> keys;
    for (u_int16_t i = 1; i < 10; i++) {
        keys.push_back("target_LA"+ std::to_string(i) + "_GHz");
    }

    for (std::string& key : keys) {
        if (parsed_input.find(key) != parsed_input.end()) {
            frequencies.push_back(parsed_input[key]);
        }
    }

    return frequencies;
}