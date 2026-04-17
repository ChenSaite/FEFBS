#include "math/fourier.h"
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <stdexcept>
#include <iostream>

namespace py = pybind11;

namespace lbcrypto {

struct __attribute__((visibility("hidden"))) FourierCalculator::Impl {
    py::object calc_func_;

    Impl() {
        py::initialize_interpreter();
        try {
            py::module_ sys = py::module_::import("sys");
            py::list path_list = sys.attr("path");
            sys.attr("path").attr("append")(PROJECT_SOURCE_DIR_STR);
            py::module_ calculator = py::module_::import("fourier_calculator");
            calc_func_             = calculator.attr("calculate_fourier_coefficients");
        }
        catch (py::error_already_set& e) {
            py::finalize_interpreter();
            throw std::runtime_error("Failed to init Python in Impl: " + std::string(e.what()));
        }
    }

    ~Impl() {
        calc_func_.release();
        py::finalize_interpreter();
    }
};

FourierCalculator::FourierCalculator() : pimpl(std::make_unique<Impl>()) {}

FourierCalculator::~FourierCalculator() = default;

FourierCalculator::FourierCalculator(FourierCalculator&&) noexcept            = default;
FourierCalculator& FourierCalculator::operator=(FourierCalculator&&) noexcept = default;

std::vector<std::complex<double>> FourierCalculator::calculate(const std::string& func_str, double a, double b,
                                                               uint32_t degree, uint32_t speed) {
    try {
        py::object result = pimpl->calc_func_(func_str, a, b, degree, speed);
        return result.cast<std::vector<std::complex<double>>>();
    }
    catch (py::error_already_set& e) {
        throw std::runtime_error("Error during Python calculation: " + std::string(e.what()));
    }
}

}  // namespace lbcrypto