# Sparse LT Integration Notes

This branch extends the FEFBS OpenFHE fork with a compiled sparse linear
transform path based on AlexanderViand's sparse LT prototype and the upstream
OpenFHE bootstrapping LT schedule.

## Implemented

- Added a resident `CompiledSparseLinearTransform` plan for sparse CKKS LT.
- Added `CompileSparseLinearTransform`,
  `EvalLinearTransformSparseCompiled`, and plan-owned rotation-key discovery.
- Reused the boot-LT double-hoist structure for baby rotations.
- Cached plaintext elements in the compiled plan to avoid repeated extraction.
- Added in-place extended ciphertext/plaintext multiply-add for group
  accumulation.
- Added guarded Horner giant accumulation when the stride key set is smaller.
- Folded direct giant first-element correction through the fast-rotation
  `addFirst` path.
- Changed extended ciphertext/plaintext multiplication initialization to use
  `CloneEmpty()` plus product elements instead of cloning full input elements.
- Added sparse LT unit tests covering direct, Horner, auto baby-step, and
  legacy evaluator parity.
- Added fixed sparse LT benchmark fixtures for AlexNet and small-model
  regression patterns.

## Deliberately not included

- No generalized fused graph compiler.
- No bootstrapping FFT-specific correction folding.
- No resident automorphism-map cache; it was not stable enough for this patch.
- No Orion-specific model or backend policy changes in this OpenFHE patch.

## Validation

Run the sparse LT contract tests:

```bash
/tmp/openfhe-sparse-lt-redesign-build/unittest/pke_tests --gtest_filter='*SparseLinearTransformCompiled*'
```

Run the fixed sparse LT benchmark set:

```bash
scripts/run_sparse_lt_benchmarks.sh 5
```

Final benchmark snapshot from this branch:

| Case | Baseline ms | Compiled ms | Ratio | Speedup | baby | keys | strategy |
|---|---:|---:|---:|---:|---:|---:|---|
| AlexNet `features_3-125` | 2577.361 | 824.181 | 0.320 | 3.13x | 2048 | 43 | direct |
| AlexNet `features_7-205` | 3021.334 | 707.307 | 0.234 | 4.27x | 256 | 58 | direct |
| LeNet `fc2-521` | 4500.847 | 1446.958 | 0.321 | 3.11x | 64 | 71 | direct |
| LoLA `conv1-13` | 444.444 | 111.012 | 0.250 | 4.00x | 2048 | 12 | direct |
| LoLA `fc2-109` | 1506.990 | 533.449 | 0.354 | 2.82x | 32 | 35 | direct |
