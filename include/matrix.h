/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <limits>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <string>
#include <cassert>

#include "vector.h"
#include "gnuplot.h"

// Data is stored by row-major so that 
// for i for j Matrix(i, j)
// iterates fast.
template<class T>
class Matrix2D {
    private:
    std::vector<std::vector<T>> data_;

    void expect_allocated() const {
        if (data_.empty()) throw std::runtime_error("Matrix is not allocated yet!");
        
        for (auto& row : data_) {
            if (!row.empty()) return;
        }

        throw std::runtime_error("Matrix is not allocated yet!");
    };

    void expect_not_allocated() const {
        if (!data_.empty()) throw std::runtime_error("Matrix is already allocated!");
    };

    public:
    Matrix2D() {};
    ~Matrix2D() {};

    // fill the vector by zeros
    void zeros(const size_t n0, const size_t n1) {
        expect_not_allocated();

        data_.resize(n0);
        for (auto& subspace : data_) {
            subspace.resize(n1); // newly added elements are value-initialized
        }
    }

    T& operator()(const size_t i, const size_t j) {
        expect_allocated();
        return data_[i][j];
    };

    const T& operator()(const size_t i, const size_t j) const {
        expect_allocated();
        return data_[i][j];
    };

    size_t nRow() const {
        expect_allocated();
        return data_.size();
    };

    size_t nCol() const {
        expect_allocated();
        return data_[0].size();
    };

    std::tuple<size_t, size_t, T> min() const {
        T value = std::numeric_limits<T>::max();
        size_t irow_min = -1;
        size_t icol_min = -1;

        for (size_t irow = 0; irow < nRow(); irow++) {
            for (size_t icol = 0; icol < nCol(); icol++) {
                const T element = operator()(irow, icol);
                if (element < value) {
                    value = element;
                    irow_min = irow;
                    icol_min = icol;
                }
            }
        }
        return std::tuple<size_t, size_t, T>(irow_min, icol_min, value);
    };

    void output_spectrum(
        const Vector<double>& x_values, 
        const Vector<double>& y_values,
        const std::string& data_file
    ) const {
        std::cout << "\nWriting spectrum data to file..." << std::endl;

        std::ofstream file(data_file);
        if (!file.is_open()) throw std::runtime_error("Could not open output file!");

        assert(x_values.size() == nRow() && "Number of data points along x-axis does not match");
        assert(y_values.size() == nCol() && "Number of data points along y-axis does not match");
        
        for (size_t i_x = 0; i_x < x_values.size(); i_x++) {
            const double x = x_values(i_x);
            for (size_t i_y = 0; i_y < y_values.size(); i_y++) {
                const double y = y_values(i_y);
                file << std::setw(7) << x;
                file << std::setw(7) << y;
                file << std::setw(14) << operator()(i_x, i_y);
                file << std::endl;
            }
            file << std::endl; // gnuplot requires a blank line after every change in the x coordinate

        }
        file.close();
    };

    // Write the data to an ASCII file (gnuplot "matrix" format) and write a
    // gnuplot script that plots it as a heat map, using user-supplied coordinates:
    //   x : coordinate of each column, x.size() must equal nCol()
    //   y : coordinate of each row,    y.size() must equal nRow()
    // The data file uses gnuplot's "matrix nonuniform" layout:
    //           x_0       x_1      ... x_{nRow-1}
    //   y_0     M[0][0]   M[1][0]  ... M[nRow-1][0]
    //   y_1     M[0][1]   M[1][1]  ... M[nRow-1][1]
    // Note that the documentation at http://gnuplot.info/docs_6.1/loc18442.html 
    // contains wrong statements:
    // "The first row of input data contains the y coordinates."
    // "The first column of input data contains the x coordinates."
    // (the first entry of the first line is ignored by gnuplot; nCol is written
    // as a placeholder).
    void write_gnuplot_xy(const Vector<T>& x,
                          const Vector<T>& y,
                          const std::string& x_label,
                          const std::string& y_label,
                          const std::string& cb_label,
                          const std::string& data_file,
                          const std::string& script_file,
                          const std::string& image_file) const {
        expect_allocated();
        if (x.size() != nCol()) throw std::invalid_argument("x.size() must be equal to nCol()");
        if (y.size() != nRow()) throw std::invalid_argument("y.size() must be equal to nRow()");

        phonon_pumping::write_heatmap_xy(data_, x.data(), y.data(), x_label, y_label, cb_label, data_file, script_file, image_file);
    }
};

// data is stored by row-major
template<class T>
class Matrix3D {
    private:
    std::vector<std::vector<std::vector<T>>> data_;

    void expect_allocated() const {
        if (data_.size() == 0u) {
            throw std::runtime_error("Matrix is not allocated yet!");
        }
    };

    void expect_not_allocated() const {
        if (data_.size() != 0u) {
            throw std::runtime_error("Matrix is already allocated!");
        }
    };

    public:
    Matrix3D() {};
    ~Matrix3D() {};

    // fill the vector by zeros
    void zeros(const size_t n0, const size_t n1, const size_t n2) {
        expect_not_allocated();

        data_.resize(n0);
        for (auto& subspace : data_) {
            subspace.resize(n1);
            for (auto& subsubspace : subspace) {
                subsubspace.resize(n2); // newly added elements are value-initialized
            }
        }
    }

    T& operator()(const size_t i, const size_t j, const size_t k) {
        expect_allocated();
        return data_[i][j][k];
    };

    const T& operator()(const size_t i, const size_t j, const size_t k) const {
        expect_allocated();
        return data_[i][j][k];
    };
};