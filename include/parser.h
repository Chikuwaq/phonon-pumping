#pragma once

#include <string>
#include <filesystem>
#include <fstream>
#include <map>

namespace phonon_pumping {

    std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        size_t end = s.find_last_not_of(" \t\r\n");
        return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    }

    std::map<std::string, double> parse_file(const std::filesystem::path& filepath) {
        if (!std::filesystem::exists(filepath)) {
            const std::string msg = "File " + filepath.string() + " not found!";
            throw std::runtime_error(msg);
        }

        std::ifstream infile(filepath);
        std::map<std::string, double> result;
        std::string line;
    
        while (std::getline(infile, line)) {
            // Remove leading/trailing whitespace
            line = trim(line);
    
            // Skip empty lines or full-line comments
            if (line.empty() || line[0] == '#')
                continue;
    
            // Remove inline comments
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) {
                line = line.substr(0, comment_pos);
                line = trim(line);
            }
    
            // Find the '=' separator
            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos)
                continue; // No key=value pattern
    
            // Extract and trim key and value
            std::string key = trim(line.substr(0, eq_pos));
            std::string value = trim(line.substr(eq_pos + 1));
    
            if (key.empty()) {
                throw std::runtime_error("Key on the left of '=' sign is empty!");
            }

            result[key] = stod(value);
        }
    
        return result;
    }

}