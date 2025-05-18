#include <functional>

namespace numerics {

    double secant_method(std::function<double(double)> f, double x0, double x1, double tol = 1e-6, int max_iter = 100);
}