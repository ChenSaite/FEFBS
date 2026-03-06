#define PROFILE

#include "openfhe.h"
#include <functional>
#include "fbs_utils.h"
using namespace lbcrypto;

int main(int argc, char* argv[]) {
    std::string boot = "x";
    std::function<double(double)> bootfunc = [](double x) -> double {
        return x;
    };

    std::string filename = "";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-f" || arg == "--file") && i + 1 < argc) {
            filename = argv[i + 1];
            i++;
        }
    }

    FuncBootstrapExample(boot, -0.5, 0.5, 25, 18, bootfunc, filename);
    return 0;
}