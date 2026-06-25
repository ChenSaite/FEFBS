#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${OPENFHE_BUILD_DIR:-/tmp/openfhe-sparse-lt-redesign-build}"
ITERATIONS="${1:-5}"
BENCH="${BUILD_DIR}/bin/examples/pke/sparse_lt_compiled_benchmark"
CASE_DIR="${SRC_DIR}/benchmark/sparse_lt_cases"
RESULT_DIR="${SRC_DIR}/benchmark/sparse_lt_results"
RESULT_FILE="${RESULT_DIR}/sparse_lt_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "${RESULT_DIR}"

cmake --build "${BUILD_DIR}" --target sparse_lt_compiled_benchmark -j "${OPENFHE_BUILD_JOBS:-2}"

{
    echo "# sparse LT benchmark"
    echo "# source=${SRC_DIR}"
    echo "# build=${BUILD_DIR}"
    echo "# iterations=${ITERATIONS}"
    echo "# timestamp=$(date -Iseconds)"

    echo "## alexnet_features_3_125_auto"
    "${BENCH}" "${CASE_DIR}/alexnet_features_3_125_auto.tsv" "${ITERATIONS}"

    echo "## alexnet_features_7_205_auto"
    "${BENCH}" "${CASE_DIR}/alexnet_features_7_205_auto.tsv" "${ITERATIONS}"

    echo "## small_regression_auto"
    "${BENCH}" "${CASE_DIR}/small_regression_auto.tsv" "${ITERATIONS}"
} | tee "${RESULT_FILE}"

echo "wrote ${RESULT_FILE}"
