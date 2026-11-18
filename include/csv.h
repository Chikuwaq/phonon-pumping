#pragma once

#include <string>
#include <filesystem>

#include "matrix.h"

namespace phonon_pumping {

    std::vector<std::string> split_csv_line(const std::string& line);
    bool read_csv(const std::filesystem::path& path, Matrix2D<std::string>& datatable);
    std::vector<double> extract_column(const Matrix2D<std::string>& datatable, const std::string& column_title);
    
}