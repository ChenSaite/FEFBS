#define PROFILE

#include "openfhe.h"
#include <functional>
#include "fbs_utils.h"
using namespace lbcrypto;

int main(int argc, char* argv[]) {
    std::string exp = "exp(x)";
    std::function<double(double)> expfunc = [](double x) -> double {
        return std::exp(x);
    };
    
    std::string filename = "";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-f" || arg == "--file") && i + 1 < argc) {
            filename = argv[i + 1];
            i++;
        }
    }
    FuncBootstrapExample(exp, -2, 2, 29, 18, expfunc, filename);
    return 0;
}