#include <fstream>
#include <iostream>
#include <sstream>

#include "csv.h"


using namespace std;

// Splits one CSV line into fields; handles quotes, escaped quotes ("") and
// keeps empty trailing fields.
std::vector<std::string> phonon_pumping::split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field += '"';   // escaped quote
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                field += c;
            }
        } else if (c == '"') {
            in_quotes = true;
        } else if (c == ',') {
            fields.push_back(std::move(field));
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(std::move(field));  // last field, even if empty
    return fields;
}

bool phonon_pumping::read_csv(
    const filesystem::path& path, 
    Matrix2D<string>& datatable
) {
    ifstream file(path); // closes the file in its destructor, so it is closed automatically when file goes out of scope at the end of the function. That also covers the early return false and any exception thrown while reading
    if (!file.is_open()) throw runtime_error("Could not open CSV file!");

    std::vector<std::vector<std::string>> rows;
    std::size_t max_cols = 0;
    std::string line;
    bool first_line = true;

    while (std::getline(file, line)) {
        // Strip UTF-8 BOM on the first line
        if (first_line) {
            if (line.size() >= 3 &&
                static_cast<unsigned char>(line[0]) == 0xEF &&
                static_cast<unsigned char>(line[1]) == 0xBB &&
                static_cast<unsigned char>(line[2]) == 0xBF) {
                line.erase(0, 3);
            }
            first_line = false;
        }

        // Strip CR from CRLF line endings
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Skip blank lines
        if (line.empty()) {
            continue;
        }

        rows.push_back(split_csv_line(line));
        max_cols = std::max(max_cols, rows.back().size());
    }

    if (file.bad()) {
        throw std::runtime_error("I/O error while reading CSV file: " + path.string());
    }

    if (rows.empty()) {
        return false;   // file opened but contained no data
    }

    datatable.zeros(rows.size(), max_cols);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (std::size_t c = 0; c < rows[r].size(); ++c) {
            datatable(r, c) = std::move(rows[r][c]);
        }
        // Short rows keep empty-string cells
    }

    return true;
}

vector<double> phonon_pumping::extract_column(
    const Matrix2D<std::string>& datatable,
    const std::string& column_title
) {
    vector<double> result;
    result.reserve(datatable.nRow());

    for (uint16_t i_col = 0; i_col < datatable.nCol(); i_col++) {
        if (datatable(0, i_col) != column_title) continue;

        for (uint16_t i_row = 0; i_row < datatable.nRow(); i_row++) {
            result.push_back(atof(datatable(i_row, i_col).c_str()));
        }
    }
    return result;
}