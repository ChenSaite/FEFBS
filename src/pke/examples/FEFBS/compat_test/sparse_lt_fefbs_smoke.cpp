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

int main(int argc, char* argv[]) {
    std::string coeffFile = "../coeffs/bootcoeff.txt";
    if (argc > 1) {
        coeffFile = argv[1];
    }

    constexpr uint32_t logn = 16;
    constexpr uint32_t slots = 1024;
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

    std::vector<std::vector<std::complex<double>>> diagonals(
        1, std::vector<std::complex<double>>(slots, std::complex<double>(1.0, 0.0)));
    std::vector<int32_t> diagonalIndices{0};
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    const uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
    const uint32_t inputQTowers = ctxt->GetElements()[0].GetNumOfElements();
    const uint32_t inputL = (inputQTowers > compositeDegree) ? inputQTowers - compositeDegree : 0;
    std::cerr << "input GetLevel=" << ctxt->GetLevel()
              << " q_towers=" << inputQTowers
              << " composite_degree=" << compositeDegree
              << " precompute_L=" << inputL << std::endl;
    auto sparseIdentity = scheme->EvalLinearTransformPrecomputeSparse(
        *cc, diagonals, diagonalIndices, 1, 1.0, inputL);

    ConstCiphertext<DCRTPoly> constInput = ctxt;
    auto afterLt = scheme->EvalLinearTransformSparse(sparseIdentity, constInput, diagonalIndices, 1);

    auto coeffs = LoadCoeffs(coeffFile);
    auto afterFefbs = cc->EvalFEFuncBootstrap(afterLt, coeffs);

    const uint32_t afterFefbsQTowers = afterFefbs->GetElements()[0].GetNumOfElements();
    const uint32_t afterFefbsL =
        (afterFefbsQTowers > compositeDegree)
            ? afterFefbsQTowers - compositeDegree
            : 0;
    std::cerr << "afterFefbs GetLevel=" << afterFefbs->GetLevel()
              << " q_towers=" << afterFefbsQTowers
              << " composite_degree=" << compositeDegree
              << " precompute_L=" << afterFefbsL << std::endl;
    auto sparseIdentityAfterFefbs = scheme->EvalLinearTransformPrecomputeSparse(
        *cc, diagonals, diagonalIndices, 1, 1.0, afterFefbsL);

    ConstCiphertext<DCRTPoly> constAfterFefbs = afterFefbs;
    auto afterSecondLt = scheme->EvalLinearTransformSparse(sparseIdentityAfterFefbs, constAfterFefbs, diagonalIndices, 1);

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, afterSecondLt, &result);
    result->SetLength(slots);
    auto got = result->GetRealPackedValue();

    double maxAbs = 0.0;
    double meanAbs = 0.0;
    for (uint32_t i = 0; i < slots; ++i) {
        double err = std::abs(got[i] - x[i]);
        maxAbs = std::max(maxAbs, err);
        meanAbs += err;
    }
    meanAbs /= static_cast<double>(slots);

    std::cout << "Sparse LT -> FEFBS identity -> Sparse LT smoke" << std::endl;
    std::cout << "level_before=" << (depth - ctxt->GetLevel()) << std::endl;
    std::cout << "level_after=" << (depth - afterSecondLt->GetLevel()) << std::endl;
    std::cout << "max_abs_error=" << maxAbs << std::endl;
    std::cout << "mean_abs_error=" << meanAbs << std::endl;
    std::cout << "precision=" << result->GetOutputPrecision(x) << " bits" << std::endl;

    if (maxAbs > 1e-4) {
        std::cerr << "combined sparse LT/FEFBS smoke failed" << std::endl;
        return 2;
    }
    return 0;
}
