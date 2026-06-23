#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-scheme.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

using namespace lbcrypto;

int main() {
    constexpr uint32_t logn = 13;
    constexpr uint32_t slots = 64;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetCKKSDataType(COMPLEX);
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << logn);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(45);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetFirstModSize(60);
    parameters.SetBatchSize(slots);
    parameters.SetMultiplicativeDepth(8);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int32_t> rotationIndices;
    for (uint32_t i = 1; i < slots; ++i) {
        rotationIndices.push_back(static_cast<int32_t>(i));
    }
    cc->EvalRotateKeyGen(keyPair.secretKey, rotationIndices);

    std::vector<double> x(slots);
    for (uint32_t i = 0; i < slots; ++i) {
        x[i] = static_cast<double>(i) / 17.0 - 1.0;
    }

    auto ptxt = cc->MakeCKKSPackedPlaintext(x, 1, 0, nullptr, slots);
    auto ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    if (!scheme) {
        std::cerr << "failed to access CKKSRNS scheme" << std::endl;
        return 1;
    }

    std::vector<int32_t> diagonalIndices{0, 1, -3, 17};
    std::vector<double> coefficients{1.0, 2.0, -0.5, 0.25};
    std::vector<std::vector<std::complex<double>>> diagonals;
    diagonals.reserve(coefficients.size());
    for (double coefficient : coefficients) {
        diagonals.emplace_back(slots, std::complex<double>(coefficient, 0.0));
    }

    auto precomputed = scheme->EvalLinearTransformPrecomputeSparse(
        *cc, diagonals, diagonalIndices, 8, 1.0, 0);
    ConstCiphertext<DCRTPoly> constCtxt = ctxt;
    auto resultCtxt = scheme->EvalLinearTransformSparse(precomputed, constCtxt, diagonalIndices, 8);

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, resultCtxt, &result);
    result->SetLength(slots);
    auto got = result->GetRealPackedValue();

    std::vector<double> expectedPlus(slots, 0.0);
    std::vector<double> expectedMinus(slots, 0.0);
    for (uint32_t i = 0; i < slots; ++i) {
        for (size_t j = 0; j < diagonalIndices.size(); ++j) {
            auto plusIndex = static_cast<uint32_t>((static_cast<int32_t>(i) + diagonalIndices[j] + slots) % slots);
            auto minusIndex = static_cast<uint32_t>((static_cast<int32_t>(i) - diagonalIndices[j] + slots) % slots);
            expectedPlus[i] += coefficients[j] * x[plusIndex];
            expectedMinus[i] += coefficients[j] * x[minusIndex];
        }
    }

    double maxAbsPlus = 0.0;
    double meanAbsPlus = 0.0;
    double maxAbsMinus = 0.0;
    double meanAbsMinus = 0.0;
    for (uint32_t i = 0; i < slots; ++i) {
        double errPlus = std::abs(got[i] - expectedPlus[i]);
        double errMinus = std::abs(got[i] - expectedMinus[i]);
        maxAbsPlus = std::max(maxAbsPlus, errPlus);
        meanAbsPlus += errPlus;
        maxAbsMinus = std::max(maxAbsMinus, errMinus);
        meanAbsMinus += errMinus;
    }
    meanAbsPlus /= static_cast<double>(slots);
    meanAbsMinus /= static_cast<double>(slots);
    const bool plusConvention = maxAbsPlus <= maxAbsMinus;
    const auto& expected = plusConvention ? expectedPlus : expectedMinus;
    const double maxAbs = plusConvention ? maxAbsPlus : maxAbsMinus;
    const double meanAbs = plusConvention ? meanAbsPlus : meanAbsMinus;

    std::cout << "Sparse LT non-identity smoke" << std::endl;
    std::cout << "diagonal_convention=" << (plusConvention ? "x[i+d]" : "x[i-d]") << std::endl;
    std::cout << "max_abs_error=" << maxAbs << std::endl;
    std::cout << "mean_abs_error=" << meanAbs << std::endl;
    std::cout << "precision=" << result->GetOutputPrecision(expected) << " bits" << std::endl;

    if (maxAbs > 1e-5) {
        std::cerr << "sparse LT non-identity smoke failed" << std::endl;
        return 2;
    }
    return 0;
}
