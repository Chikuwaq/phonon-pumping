#pragma once

#include <cstdio>
#include <vector>
#include <stdexcept>

template<class T>
class Vector {
    private:
    std::vector<T> data_;

    public:
    Vector() {};
    ~Vector() {};

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

        data_.resize(nElements);
        for (T& element : data_) {
            element = static_cast<T>(0);
        }
    }

    T& operator()(const size_t i) {
        return data_[i];
    };

    const T& operator()(const size_t i) const {
        return data_[i];
    };

    size_t size() const {
        return data_.size();
    }
};