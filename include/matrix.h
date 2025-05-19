/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdio>
#include <vector>

// Data is stored by row-major so that 
// for i for j Matrix(i, j)
// iterates fast.
template<class T>
class Matrix2D {
    private:
    std::vector<std::vector<T>> data_;

    void check_allocated() const {
        if (data_.size() == 0u) {
            throw std::runtime_error("Matrix is not allocated yet!");
        }
    };

    void check_not_allocated() const {
        if (data_.size() != 0u) {
            throw std::runtime_error("Matrix is already allocated!");
        }
    };

    public:
    Matrix2D() {};
    ~Matrix2D() {};

    // fill the vector by zeros
    void zeros(const size_t n0, const size_t n1) {
        check_not_allocated();

        data_.resize(n0);
        for (auto& subspace : data_) {
            subspace.resize(n1);
            for (T& element : subspace) {
                element = static_cast<T>(0);
            }
        }
    }

    T& operator()(const size_t i, const size_t j) {
        check_allocated();
        return data_[i][j];
    };

    const T& operator()(const size_t i, const size_t j) const {
        check_allocated();
        return data_[i][j];
    };

    size_t nRow() const {
        check_allocated();
        return data_.size();
    };

    size_t nCol() const {
        check_allocated();
        return data_[0].size();
    };
};

// data is stored by row-major
template<class T>
class Matrix3D {
    private:
    std::vector<std::vector<std::vector<T>>> data_;

    void check_allocated() const {
        if (data_.size() == 0u) {
            throw std::runtime_error("Matrix is not allocated yet!");
        }
    };

    void check_not_allocated() const {
        if (data_.size() != 0u) {
            throw std::runtime_error("Matrix is already allocated!");
        }
    };

    public:
    Matrix3D() {};
    ~Matrix3D() {};

    // fill the vector by zeros
    void zeros(const size_t n0, const size_t n1, const size_t n2) {
        check_not_allocated();

        data_.resize(n0);
        for (auto& subspace : data_) {
            subspace.resize(n1);
            for (auto& subsubspace : subspace) {
                subsubspace.resize(n2);
                for (T& element : subsubspace) {
                    element = static_cast<T>(0);
                }
            }
        }
    }

    T& operator()(const size_t i, const size_t j, const size_t k) {
        check_allocated();
        return data_[i][j][k];
    };

    const T& operator()(const size_t i, const size_t j, const size_t k) const {
        check_allocated();
        return data_[i][j][k];
    };
};