// Contract test for the resident sparse-LT execution plan.  It compares the
// legacy sparse evaluator with the bootstrap-style compiled schedule using
// identical pre-encoded plaintext diagonals.

#include "gtest/gtest.h"
#include "openfhe.h"

#include <algorithm>
#include <complex>
#include <set>
#include <vector>

using namespace lbcrypto;

class UnitTestSparseLinearTransformCompiled : public ::testing::TestWithParam<std::vector<int32_t>> {};

static void ExpectNearPlaintexts(const std::vector<std::complex<double>>& expected,
                                 const std::vector<std::complex<double>>& actual,
                                 double tolerance) {
    ASSERT_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_NEAR(expected[i].real(), actual[i].real(), tolerance) << "slot " << i;
}

static std::set<int32_t> ExpectedCompiledRotations(const CompiledSparseLinearTransform& plan) {
    std::set<int32_t> rotations;
    for (const auto baby : plan.babyRotations)
        rotations.insert(static_cast<int32_t>(baby));
    if (plan.useHornerGiantAccumulation) {
        for (const auto giant : plan.hornerGiantRotations)
            rotations.insert(static_cast<int32_t>(giant));
    }
    else {
        for (const auto& group : plan.groups) {
            if (group.giantRotation != 0)
                rotations.insert(static_cast<int32_t>(group.giantRotation));
        }
    }
    return rotations;
}

TEST_P(UnitTestSparseLinearTransformCompiled, MatchesLegacySparseEvaluator) {
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(4);
    parameters.SetScalingModSize(50);
    parameters.SetBatchSize(8);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    // Sparse LT reduces negative diagonals to their positive slot rotations,
    // so install the complete BSGS rotation domain for this contract test.
    cc->EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4, 5, 6, 7});

    const std::vector<std::complex<double>> input{{0.125, 0.0}, {0.25, 0.0}, {0.375, 0.0}, {0.5, 0.0},
                                                   {0.625, 0.0}, {0.75, 0.0}, {0.875, 0.0}, {1.0, 0.0}};
    auto ciphertext = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(input));

    const auto& indices = GetParam();
    std::vector<std::vector<std::complex<double>>> diagonals(
        indices.size(), std::vector<std::complex<double>>(input.size(), {0.0, 0.0}));
    for (size_t d = 0; d < diagonals.size(); ++d)
        for (size_t s = 0; s < input.size(); ++s)
            diagonals[d][s] = {0.0625 * static_cast<double>(d + 1), 0.0};

    constexpr uint32_t kBabyStep = 2;
    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    ASSERT_NE(scheme, nullptr);
    auto precomputed = scheme->EvalLinearTransformPrecomputeSparse(*cc, diagonals, indices, kBabyStep);
    ConstCiphertext<DCRTPoly> constCiphertext = ciphertext;
    auto legacy = scheme->EvalLinearTransformSparse(precomputed, constCiphertext, indices, kBabyStep);
    auto plan = scheme->CompileSparseLinearTransform(precomputed, indices, ciphertext->GetSlots(), kBabyStep,
                                                     cc->GetCyclotomicOrder());
    auto compiled = scheme->EvalLinearTransformSparseCompiled(*plan, constCiphertext);
    const auto& planRotations = scheme->GetSparseLinearTransformRotationIndices(*plan);
    EXPECT_EQ(std::set<int32_t>(planRotations.begin(), planRotations.end()), ExpectedCompiledRotations(*plan));

    Plaintext legacyPlain;
    Plaintext compiledPlain;
    cc->Decrypt(keys.secretKey, legacy, &legacyPlain);
    cc->Decrypt(keys.secretKey, compiled, &compiledPlain);
    legacyPlain->SetLength(input.size());
    compiledPlain->SetLength(input.size());
    ExpectNearPlaintexts(legacyPlain->GetCKKSPackedValue(), compiledPlain->GetCKKSPackedValue(), 1e-6);
}

INSTANTIATE_TEST_SUITE_P(SparseSchedules, UnitTestSparseLinearTransformCompiled,
                         ::testing::Values(std::vector<int32_t>{0, 1, 2, 3},
                                           std::vector<int32_t>{-3, -1, 0, 2},
                                           std::vector<int32_t>{0, 2},
                                           std::vector<int32_t>{0, 1}));

TEST(SparseLinearTransformCompiled, AutoBabyStepUsesSparseWorkSchedule) {
    constexpr uint32_t kSlots = 16;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(4);
    parameters.SetScalingModSize(50);
    parameters.SetBatchSize(kSlots);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    std::vector<int32_t> rotations(kSlots - 1);
    for (uint32_t i = 1; i < kSlots; ++i)
        rotations[i - 1] = static_cast<int32_t>(i);
    cc->EvalRotateKeyGen(keys.secretKey, rotations);

    std::vector<std::complex<double>> input(kSlots);
    for (uint32_t i = 0; i < kSlots; ++i)
        input[i] = {static_cast<double>(i + 1) / 32.0, 0.0};
    auto ciphertext = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(input));

    const std::vector<int32_t> indices{0, 1, 2, 3, 4};
    std::vector<std::vector<std::complex<double>>> diagonals(indices.size(),
                                                             std::vector<std::complex<double>>(kSlots));
    for (size_t d = 0; d < diagonals.size(); ++d)
        for (uint32_t slot = 0; slot < kSlots; ++slot)
            diagonals[d][slot] = {0.015625 * static_cast<double>((d + 1) * ((slot % 2) + 1)), 0.0};

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    ASSERT_NE(scheme, nullptr);
    auto plaintexts = scheme->EvalLinearTransformPrecomputeSparse(*cc, diagonals, indices, 0);
    auto plan = scheme->CompileSparseLinearTransform(plaintexts, indices, kSlots, 0, cc->GetCyclotomicOrder());
    EXPECT_GT(plan->babyStep, 0u);
    EXPECT_LE(plan->babyStep, kSlots);
    EXPECT_EQ(plan->babyStep & (plan->babyStep - 1), 0u);
    const auto& planRotations = scheme->GetSparseLinearTransformRotationIndices(*plan);
    EXPECT_EQ(std::set<int32_t>(planRotations.begin(), planRotations.end()), ExpectedCompiledRotations(*plan));

    ConstCiphertext<DCRTPoly> inputCt = ciphertext;
    auto legacy = scheme->EvalLinearTransformSparse(plaintexts, inputCt, indices, 0);
    auto compiled = scheme->EvalLinearTransformSparseCompiled(*plan, inputCt);

    Plaintext legacyPlain;
    Plaintext compiledPlain;
    cc->Decrypt(keys.secretKey, legacy, &legacyPlain);
    cc->Decrypt(keys.secretKey, compiled, &compiledPlain);
    legacyPlain->SetLength(kSlots);
    compiledPlain->SetLength(kSlots);
    ExpectNearPlaintexts(legacyPlain->GetCKKSPackedValue(), compiledPlain->GetCKKSPackedValue(), 1e-6);
}

TEST(SparseLinearTransformCompiled, UsesHornerGiantAccumulationWhenStrideKeySetIsSmaller) {
    constexpr uint32_t kSlots = 64;
    constexpr uint32_t kBabyStep = 4;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(4);
    parameters.SetScalingModSize(50);
    parameters.SetBatchSize(kSlots);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    std::vector<int32_t> rotations(kSlots - 1);
    for (uint32_t i = 1; i < kSlots; ++i)
        rotations[i - 1] = static_cast<int32_t>(i);
    cc->EvalRotateKeyGen(keys.secretKey, rotations);

    std::vector<std::complex<double>> input(kSlots);
    for (uint32_t i = 0; i < kSlots; ++i)
        input[i] = {static_cast<double>((i % 11) + 1) / 128.0, 0.0};
    auto ciphertext = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(input));

    std::vector<int32_t> indices;
    for (uint32_t giant = 0; giant <= 44; giant += kBabyStep)
        indices.push_back(static_cast<int32_t>(giant));

    std::vector<std::vector<std::complex<double>>> diagonals(indices.size(),
                                                             std::vector<std::complex<double>>(kSlots));
    for (size_t d = 0; d < diagonals.size(); ++d)
        for (uint32_t slot = 0; slot < kSlots; ++slot)
            diagonals[d][slot] = {0.00390625 * static_cast<double>((d + 1) * ((slot % 5) + 1)), 0.0};

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    ASSERT_NE(scheme, nullptr);
    auto plaintexts = scheme->EvalLinearTransformPrecomputeSparse(*cc, diagonals, indices, kBabyStep);
    ConstCiphertext<DCRTPoly> inputCt = ciphertext;
    auto legacy = scheme->EvalLinearTransformSparse(plaintexts, inputCt, indices, kBabyStep);
    auto plan = scheme->CompileSparseLinearTransform(plaintexts, indices, kSlots, kBabyStep, cc->GetCyclotomicOrder());
    auto compiled = scheme->EvalLinearTransformSparseCompiled(*plan, inputCt);

    EXPECT_TRUE(plan->useHornerGiantAccumulation);
    EXPECT_EQ(plan->hornerGiantRotations, std::vector<uint32_t>{kBabyStep});
    const auto& planRotations = scheme->GetSparseLinearTransformRotationIndices(*plan);
    EXPECT_EQ(std::set<int32_t>(planRotations.begin(), planRotations.end()), ExpectedCompiledRotations(*plan));

    Plaintext legacyPlain;
    Plaintext compiledPlain;
    cc->Decrypt(keys.secretKey, legacy, &legacyPlain);
    cc->Decrypt(keys.secretKey, compiled, &compiledPlain);
    legacyPlain->SetLength(kSlots);
    compiledPlain->SetLength(kSlots);
    ExpectNearPlaintexts(legacyPlain->GetCKKSPackedValue(), compiledPlain->GetCKKSPackedValue(), 1e-5);
}

TEST(SparseLinearTransformCompiled, MatchesLegacySparseEvaluatorForWideSparsePattern) {
    constexpr uint32_t kSlots = 32;
    constexpr uint32_t kBabyStep = 8;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetCKKSDataType(COMPLEX);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 12);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(45);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetFirstModSize(60);
    parameters.SetBatchSize(kSlots);
    parameters.SetMultiplicativeDepth(8);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);
    std::vector<int32_t> rotations(kSlots - 1);
    for (uint32_t i = 1; i < kSlots; ++i)
        rotations[i - 1] = static_cast<int32_t>(i);
    cc->EvalRotateKeyGen(keys.secretKey, rotations);

    std::vector<std::complex<double>> input(kSlots);
    for (uint32_t i = 0; i < kSlots; ++i)
        input[i] = {static_cast<double>(i + 1) / 64.0, 0.0};
    auto ciphertext = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(input));

    const std::vector<int32_t> indices{0, 1, 2, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31, -1, -9};
    std::vector<std::vector<std::complex<double>>> diagonals(indices.size(),
                                                             std::vector<std::complex<double>>(kSlots));
    for (size_t d = 0; d < diagonals.size(); ++d)
        for (uint32_t slot = 0; slot < kSlots; ++slot)
            diagonals[d][slot] = {0.0078125 * static_cast<double>((d + 1) * ((slot % 3) + 1)), 0.0};

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    ASSERT_NE(scheme, nullptr);
    auto plaintexts = scheme->EvalLinearTransformPrecomputeSparse(*cc, diagonals, indices, kBabyStep);
    ConstCiphertext<DCRTPoly> inputCt = ciphertext;
    auto legacy = scheme->EvalLinearTransformSparse(plaintexts, inputCt, indices, kBabyStep);
    auto plan = scheme->CompileSparseLinearTransform(plaintexts, indices, kSlots, kBabyStep, cc->GetCyclotomicOrder());
    auto compiled = scheme->EvalLinearTransformSparseCompiled(*plan, inputCt);

    EXPECT_EQ(plan->babyStep, kBabyStep);
    const auto& planRotations = scheme->GetSparseLinearTransformRotationIndices(*plan);
    EXPECT_EQ(std::set<int32_t>(planRotations.begin(), planRotations.end()), ExpectedCompiledRotations(*plan));
    for (const auto& group : plan->groups) {
        if (group.giantRotation != 0)
            EXPECT_NE(group.automorphismIndex, 0u);
    }

    Plaintext legacyPlain;
    Plaintext compiledPlain;
    cc->Decrypt(keys.secretKey, legacy, &legacyPlain);
    cc->Decrypt(keys.secretKey, compiled, &compiledPlain);
    legacyPlain->SetLength(kSlots);
    compiledPlain->SetLength(kSlots);
    ExpectNearPlaintexts(legacyPlain->GetCKKSPackedValue(), compiledPlain->GetCKKSPackedValue(), 1e-5);
}
