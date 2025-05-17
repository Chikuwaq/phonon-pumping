#pragma once

#include <string>
#include <filesystem>
#include <fstream>
#include <map>

namespace phonon_pumping {

    std::string trim(const std::string& s);

    std::map<std::string, double> parse_file(const std::filesystem::path& filepath);

}