#define PROFILE

#include "openfhe.h"
#include <functional>
#include "fbs_utils.h"
using namespace lbcrypto;

int main(int argc, char* argv[]) {
    std::string gelu                       = "0.5*x*(1+tanh(sqrt(2/pi)*(x+0.044715*x**3)))";
    
    std::function<double(double)> gelufunc = [](double x) -> double {
        return 0.5 * x * (1 + std::tanh(std::sqrt(2 / M_PI) * (x + 0.044715 * std::pow(x, 3))));
    };

    std::string filename = "";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-f" || arg == "--file") && i + 1 < argc) {
            filename = argv[i + 1];
            i++;
        }
    }
    FuncBootstrapExample(gelu, -8, 8, 44, 12, gelufunc, filename);
    return 0;
}

