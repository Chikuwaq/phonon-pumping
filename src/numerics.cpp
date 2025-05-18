#include <iostream>

#include "numerics.h"

// find root of the function by iteration
double numerics::secant_method(std::function<double(double)> f, double x0, double x1, double tol, int max_iter) {
    double f0 = f(x0);
    double f1 = f(x1);

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::abs(f1 - f0) < 1e-12) throw std::runtime_error("Division by zero in secant method.");
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        if (std::abs(x2 - x1) < tol) {
            std::cout << "Secant method found the root at iteration " << iter << std::endl;
            return x2;
        }
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }
    throw std::runtime_error("Secant method did not converge.");
}