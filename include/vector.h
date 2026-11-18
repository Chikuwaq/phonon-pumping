/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>

template<class T>
class Vector {
    public:
    using const_iterator = const T*;

    private:
    std::vector<T> data_;

    public:
    Vector() {};
    ~Vector() {};

    Vector(const Vector<T>& other)
    : data_(other.data_) {
    }

    Vector(Vector<T>&& other) {
        data_ = std::move(other.data_);
    }

    Vector<T>& operator=(const Vector<T>& other) {
        data_ = other.data_;
    }

    Vector<T>& operator=(Vector<T>&& other) {
        data_ = std::move(other.data_);
    }

    const_iterator begin() const { 
		return data_.data(); 
	}
	
	const_iterator end() const { 
		return data_.data() + size(); 
	}

    // fill the vector by equidistant values in the range [min, max)
    void linspace(const T min, const T max, const size_t nElements) {
        if (data_.size() != 0u) {
            throw std::runtime_error("Vector is already allocated!");
        }

        data_.resize(nElements);
        for (size_t i = 0; i < nElements; i++) {
            data_[i] = (max - min) * i / nElements + min;
        }
    };

    // fill the vector by zeros
    void zeros(const size_t nElements) {
        if (data_.size() != 0u) {
            throw std::runtime_error("Vector is already allocated!");
        }

        data_.resize(nElements); // newly added elements are value-initialized
    }

    T& operator()(const size_t i) {
        return data_[i];
    };

    const T& operator()(const size_t i) const {
        return data_[i];
    };

    const std::vector<T>& data() const {
        return data_;
    }

    size_t size() const {
        return data_.size();
    }

    void println() const {
        for (const T& element : data_) {
            std::cout << element << std::endl;
        }
    }
};