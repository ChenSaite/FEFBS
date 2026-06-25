// Compare the baseline sparse LT evaluator with the bootstrap-style compiled
// double-hoist schedule.  The cases below are synthetic Orion/AlexNet-like
// sparse diagonal distributions; they do not depend on Orion or Python.
#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-scheme.h"

#include <algorithm>
#include <chrono>
#include <complex>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace lbcrypto;

template <typename F>
double MedianMilliseconds(F&& operation, uint32_t iterations) {
    std::vector<double> samples;
    samples.reserve(iterations);
    for (uint32_t i = 0; i < iterations; ++i) {
        const auto start = std::chrono::steady_clock::now();
        auto result = operation();
        const auto end = std::chrono::steady_clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(end - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

std::vector<int32_t> StridedSparseIndices(uint32_t slots, uint32_t count, uint32_t stride) {
    std::vector<int32_t> indices;
    indices.reserve(count);
    uint32_t value = 0;
    for (uint32_t i = 0; i < count; ++i) {
        indices.push_back(static_cast<int32_t>(value % slots));
        value += stride;
    }
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

std::vector<int32_t> RangesToIndices(const std::vector<std::pair<int32_t, int32_t>>& ranges) {
    std::vector<int32_t> indices;
    for (const auto& [start, end] : ranges)
        for (int32_t value = start; value <= end; ++value)
            indices.push_back(value);
    return indices;
}

std::vector<std::tuple<std::string, uint32_t, uint32_t, std::vector<int32_t>>> LoadCases(const std::string& path) {
    std::ifstream input(path);
    if (!input.is_open())
        OPENFHE_THROW("Unable to open sparse LT case file: " + path);
    std::vector<std::tuple<std::string, uint32_t, uint32_t, std::vector<int32_t>>> cases;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        std::stringstream ss(line);
        std::string name;
        std::string slotsText;
        std::string babyStepText;
        std::string rangesText;
        if (!std::getline(ss, name, '\t') || !std::getline(ss, slotsText, '\t') ||
            !std::getline(ss, babyStepText, '\t') || !std::getline(ss, rangesText, '\t')) {
            OPENFHE_THROW("Malformed sparse LT case line: " + line);
        }
        std::vector<std::pair<int32_t, int32_t>> ranges;
        std::stringstream rangesStream(rangesText);
        std::string rangeText;
        while (std::getline(rangesStream, rangeText, ',')) {
            auto dash = rangeText.find('-');
            if (dash == std::string::npos)
                OPENFHE_THROW("Malformed sparse LT range: " + rangeText);
            int32_t start = std::stoi(rangeText.substr(0, dash));
            int32_t end   = std::stoi(rangeText.substr(dash + 1));
            ranges.emplace_back(start, end);
        }
        cases.emplace_back(name, static_cast<uint32_t>(std::stoul(slotsText)),
                           static_cast<uint32_t>(std::stoul(babyStepText)),
                           RangesToIndices(ranges));
    }
    return cases;
}

std::vector<std::vector<std::complex<double>>> MakeDiagonals(const std::vector<int32_t>& indices,
                                                             uint32_t slots) {
    std::vector<std::vector<std::complex<double>>> diagonals(indices.size(),
                                                             std::vector<std::complex<double>>(slots));
    for (size_t d = 0; d < diagonals.size(); ++d)
        for (uint32_t slot = 0; slot < slots; ++slot)
            diagonals[d][slot] = {0.0009765625 * static_cast<double>((d % 17) + 1) *
                                      static_cast<double>((slot % 5) + 1),
                                  0.0};
    return diagonals;
}

void RunCase(const std::string& name, const CryptoContext<DCRTPoly>& cc,
             const PublicKey<DCRTPoly>& publicKey, const PrivateKey<DCRTPoly>& secretKey,
             std::shared_ptr<SchemeCKKSRNS> scheme,
             uint32_t slots, uint32_t babyStep, const std::vector<int32_t>& indices,
             uint32_t iterations) {
    std::vector<std::complex<double>> input(slots);
    for (uint32_t i = 0; i < slots; ++i)
        input[i] = {static_cast<double>(i) / 257.0 - 0.5, 0.0};
    auto ciphertext = cc->Encrypt(publicKey, cc->MakeCKKSPackedPlaintext(input));
    auto diagonals = MakeDiagonals(indices, slots);
    auto plaintexts = scheme->EvalLinearTransformPrecomputeSparse(*cc, diagonals, indices, babyStep);
    ConstCiphertext<DCRTPoly> inputCt = ciphertext;
    auto plan = scheme->CompileSparseLinearTransform(plaintexts, indices, slots, babyStep,
                                                     cc->GetCyclotomicOrder());
    std::set<int32_t> keySet;
    const auto& planKeySet = scheme->GetSparseLinearTransformRotationIndices(*plan);
    keySet.insert(planKeySet.begin(), planKeySet.end());
    for (const auto index : indices) {
        const uint32_t rotation = ReduceRotation(index, slots);
        const uint32_t baby     = rotation & (plan->babyStep - 1);
        const uint32_t giant    = ((rotation / plan->babyStep) * plan->babyStep) % slots;
        if (baby != 0)
            keySet.insert(static_cast<int32_t>(baby));
        if (giant != 0)
            keySet.insert(static_cast<int32_t>(giant));
    }
    std::vector<int32_t> benchmarkKeys(keySet.begin(), keySet.end());
    cc->EvalRotateKeyGen(secretKey, benchmarkKeys);

    scheme->EvalLinearTransformSparse(plaintexts, inputCt, indices, babyStep);
    scheme->EvalLinearTransformSparseCompiled(*plan, inputCt);
    const double baselineMs = MedianMilliseconds(
        [&] { return scheme->EvalLinearTransformSparse(plaintexts, inputCt, indices, babyStep); }, iterations);
    const double compiledMs = MedianMilliseconds(
        [&] { return scheme->EvalLinearTransformSparseCompiled(*plan, inputCt); }, iterations);

    size_t nonzeroGiantGroups = 0;
    for (const auto& group : plan->groups) {
        if (group.giantRotation != 0)
            ++nonzeroGiantGroups;
    }
    const auto& rotationIndices = scheme->GetSparseLinearTransformRotationIndices(*plan);
    std::cout << std::fixed << std::setprecision(3)
              << "case=" << name
              << " slots=" << slots
              << " diagonals=" << indices.size()
              << " requested_baby_step=" << babyStep
              << " actual_baby_step=" << plan->babyStep
              << " giant_strategy=" << (plan->useHornerGiantAccumulation ? "horner" : "direct")
              << " rotation_keys=" << rotationIndices.size()
              << " unique_baby=" << plan->babyRotations.size()
              << " horner_giant_keys=" << plan->hornerGiantRotations.size()
              << " giant_groups=" << plan->groups.size()
              << " nonzero_giant_groups=" << nonzeroGiantGroups
              << " baseline_median_ms=" << baselineMs
              << " compiled_median_ms=" << compiledMs
              << " compiled_over_baseline=" << (compiledMs / baselineMs) << '\n';
}

int main(int argc, char** argv) {
    constexpr uint32_t kSlots = 256;
    uint32_t iterations = 5;
    const std::string casePath = (argc > 1) ? argv[1] : "";
    if (argc > 2)
        iterations = static_cast<uint32_t>(std::stoul(argv[2]));

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetCKKSDataType(COMPLEX);
    parameters.SetSecretKeyDist(SPARSE_TERNARY);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 16);
    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(45);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    parameters.SetFirstModSize(60);
    parameters.SetBatchSize(32768);
    parameters.SetMultiplicativeDepth(8);
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    auto scheme = std::dynamic_pointer_cast<SchemeCKKSRNS>(cc->GetScheme());
    if (!scheme)
        return 1;
    if (!casePath.empty()) {
        for (const auto& [name, slots, babyStep, indices] : LoadCases(casePath))
            RunCase(name, cc, keys.publicKey, keys.secretKey, scheme, slots, babyStep, indices, iterations);
        return 0;
    }

    RunCase("synthetic-clustered-small", cc, keys.publicKey, keys.secretKey, scheme, kSlots, 16,
            std::vector<int32_t>{0, 1, 2, 3, 16, 17, 18, 32, 33, 64, 65, 128}, iterations);
    RunCase("synthetic-strided-medium", cc, keys.publicKey, keys.secretKey, scheme, kSlots, 16,
            StridedSparseIndices(kSlots, 64, 17), iterations);
    RunCase("synthetic-alexnet-like-wide", cc, keys.publicKey, keys.secretKey, scheme, kSlots, 32,
            StridedSparseIndices(kSlots, 160, 29), iterations);
    return 0;
}
