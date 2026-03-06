#define PROFILE

#include "openfhe.h"
#include <functional>
#include "fbs_utils.h"
using namespace lbcrypto;

int main(int argc, char* argv[]) {
    std::string sigmoid                    = "1/(1+exp(-x))";
    std::function<double(double)> sigmoidfunc = [](double x) -> double {
        return 1 / (1 + std::exp(-x));
    };

    std::string filename = "";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-f" || arg == "--file") && i + 1 < argc) {
            filename = argv[i + 1];
            i++;
        }
    }
    FuncBootstrapExample(sigmoid, -8, 8, 34, 9, sigmoidfunc, filename);
    return 0;
}