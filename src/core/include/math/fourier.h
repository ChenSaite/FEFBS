#ifndef LBCRYPTO_MATH_FOURIER_H
#define LBCRYPTO_MATH_FOURIER_H

#include <vector>
#include <complex>
#include <string>
#include <memory>

namespace lbcrypto {
class FourierCalculator {
public:
    // --- Declarations ---
    FourierCalculator();
    ~FourierCalculator(); // The destructor must be declared

    // Declare copy and move operations
    FourierCalculator(const FourierCalculator&) = delete;
    FourierCalculator& operator=(const FourierCalculator&) = delete;
    FourierCalculator(FourierCalculator&&) noexcept;
    FourierCalculator& operator=(FourierCalculator&&) noexcept;

    std::vector<std::complex<double>> calculate(
        const std::string& func_str, double a, double b,
        std::uint32_t degree, std::uint32_t speed = 1);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};
const char* PROJECT_SOURCE_DIR_STR = "@PROJECT_SOURCE_DIR@";
}

#endif // LBCRYPTO_MATH_FOURIER_H