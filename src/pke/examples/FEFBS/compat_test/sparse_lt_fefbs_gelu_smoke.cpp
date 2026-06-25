#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-scheme.h"
#include "../fbs_utils.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

using namespace lbcrypto;

namespace {
double Gelu(double x) {
    return 0.5 * x * (1.0 + std::tanh(std::sqrt(2 / M_PI) * (x + 0.044715 * std::pow(x, 3))));
}

std::vector<double> ApplyTwoDiagonalLt(const std::vector<double>& x, double c0, double c1) {
    const uint32_t slots = x.size();
    std::vector<double> y(slots);
    for (uint32_t i = 0; i < slots; ++i) {
        y[i] = c0 * x[i] + c1 * x[(i + 1) % slots];
    }
    return y;
}

uint32_t LevelForCiphertext(const CryptoContext<DCRTPoly>& cc, ConstCiphertext<DCRTPoly> ctxt) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    const uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
    const uint32_t qTowers = ctxt->GetElements()[0].GetNumOfElements();
    return (qTowers > compositeDegree) ? qTowers - compositeDegree : 0;
}
}  // namespace

int main(int argc, char* argv[]) {
    std::string coeffFile = "../coeffs/gelucoeff.txt";
    if (argc > 1) {
        coeffFile = argv[1];
    }

    constexpr uint32_t logn = 16;
    constexpr uint32_t slots = 256;
    std::vector<uint32_t> levelBudget = {3, 2};
    std::vector<uint32_t> bsgsDim = {0, 0};
    const uint32_t depth = levelBudget[0] + levelBudget[1] + 12 + 9;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetCKKSDataType(COMPLEX);
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << logn);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(59);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetFirstModSize(60);
    parameters.SetBatchSize(slots);
    parameters.SetMultiplicativeDepth(depth);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cc->EvalFEFuncBootstrapSetup(levelBudget, bsgsDim, slots);
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, slots);

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    if (!scheme) {
        std::cerr << "failed to access CKKSRNS scheme" << std::endl;
        return 1;
    }

    std::vector<double> x(slots);
    for (uint32_t i = 0; i < slots; ++i) {
        x[i] = -0.5 + static_cast<double>(i) / static_cast<double>(slots);
    }

    auto ptxt = cc->MakeCKKSPackedPlaintext(x, 1, depth - (levelBudget[1] + 1), nullptr, slots);
    auto ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    std::vector<int32_t> diagonalIndices{0, 1};
    std::vector<std::vector<std::complex<double>>> firstDiagonals{
        std::vector<std::complex<double>>(slots, std::complex<double>(0.7, 0.0)),
        std::vector<std::complex<double>>(slots, std::complex<double>(0.2, 0.0)),
    };
    std::vector<std::vector<std::complex<double>>> secondDiagonals{
        std::vector<std::complex<double>>(slots, std::complex<double>(1.1, 0.0)),
        std::vector<std::complex<double>>(slots, std::complex<double>(-0.1, 0.0)),
    };

    ConstCiphertext<DCRTPoly> constInput = ctxt;
    auto firstLtPrecompute = scheme->EvalLinearTransformPrecomputeSparse(
        *cc, firstDiagonals, diagonalIndices, 1, 1.0, LevelForCiphertext(cc, constInput));
    auto firstLtPlan = scheme->CompileSparseLinearTransform(
        firstLtPrecompute, diagonalIndices, slots, 1, cc->GetCyclotomicOrder());
    const auto& firstLtKeys = scheme->GetSparseLinearTransformRotationIndices(*firstLtPlan);
    if (!firstLtKeys.empty())
        cc->EvalRotateKeyGen(keyPair.secretKey, firstLtKeys);
    auto afterFirstLt = scheme->EvalLinearTransformSparseCompiled(*firstLtPlan, constInput);

    auto coeffs = LoadCoeffs(coeffFile);
    auto afterGeluFefbs = cc->EvalFEFuncBootstrap(afterFirstLt, coeffs);

    ConstCiphertext<DCRTPoly> constAfterGelu = afterGeluFefbs;
    auto secondLtPrecompute = scheme->EvalLinearTransformPrecomputeSparse(
        *cc, secondDiagonals, diagonalIndices, 1, 1.0, LevelForCiphertext(cc, constAfterGelu));
    auto secondLtPlan = scheme->CompileSparseLinearTransform(
        secondLtPrecompute, diagonalIndices, slots, 1, cc->GetCyclotomicOrder());
    const auto& secondLtKeys = scheme->GetSparseLinearTransformRotationIndices(*secondLtPlan);
    if (!secondLtKeys.empty())
        cc->EvalRotateKeyGen(keyPair.secretKey, secondLtKeys);
    auto resultCtxt = scheme->EvalLinearTransformSparseCompiled(*secondLtPlan, constAfterGelu);

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, resultCtxt, &result);
    result->SetLength(slots);
    auto got = result->GetRealPackedValue();

    auto firstExpected = ApplyTwoDiagonalLt(x, 0.7, 0.2);
    std::vector<double> geluExpected(slots);
    for (uint32_t i = 0; i < slots; ++i) {
        geluExpected[i] = Gelu(16.0 * firstExpected[i]);
    }
    auto expected = ApplyTwoDiagonalLt(geluExpected, 1.1, -0.1);

    double maxAbs = 0.0;
    double meanAbs = 0.0;
    for (uint32_t i = 0; i < slots; ++i) {
        double err = std::abs(got[i] - expected[i]);
        maxAbs = std::max(maxAbs, err);
        meanAbs += err;
    }
    meanAbs /= static_cast<double>(slots);

    std::cout << "Sparse LT -> FEFBS GELU -> Sparse LT smoke" << std::endl;
    std::cout << "level_before=" << (depth - ctxt->GetLevel()) << std::endl;
    std::cout << "level_after_first_lt=" << (depth - afterFirstLt->GetLevel()) << std::endl;
    std::cout << "level_after_gelu_fefbs=" << (depth - afterGeluFefbs->GetLevel()) << std::endl;
    std::cout << "level_after_second_lt=" << (depth - resultCtxt->GetLevel()) << std::endl;
    std::cout << "max_abs_error=" << maxAbs << std::endl;
    std::cout << "mean_abs_error=" << meanAbs << std::endl;
    std::cout << "precision=" << result->GetOutputPrecision(expected) << " bits" << std::endl;

    if (maxAbs > 1e-3) {
        std::cerr << "combined sparse LT/FEFBS GELU smoke failed" << std::endl;
        return 2;
    }
    return 0;
}
