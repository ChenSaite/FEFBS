# High-Precision Functional Bootstrapping for CKKS from Fourier Extension
This repository contains the implementation of the paper ["High-Precision Functional Bootstrapping for CKKS from Fourier Extension"](https://eprint.iacr.org/2026/367), built on the [OpenFHE](https://github.com/openfheorg/openfhe-development) library.

## Installation
### System Requirements
Please refer to [OpenFHE Installation Documentation](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/linux.html).
### Python Dependencies
The Fourier coefficient generation relies on Python scripts. Please ensure you have Python 3 installed with the following libraries:
```bash
pip install numpy sympy
```
### Building the Project
```bash
mkdir build
cd build
cmake ..
make -j
```

## Homomorphic Encryption Benchmarking
The source code for the test experiments is located in `src/pke/examples/FEFBS/`, which includes test files various benchmarking examples. (e.g. GeLU(x), exp(x)) 

After a successful build, the corresponding executable binaries are generated in the `build/bin/examples/pke/` directory.

The test binaries support two modes of operation for obtaining the necessary Fourier coefficients:
### 1. Online Coefficient Generation (Python Dependency)
Our framework invokes `fourier_calculator.py` at runtime to calculate Fourier coefficients. Therefore, you must specify the PYTHONPATH environment variable so the binary can locate the script.

To run a test with online generation (e.g., the exp function):
```bash
# Replace [script_path] with the directory of fourier_calculator.py
PYTHONPATH=[script_path] ./build/bin/examples/pke/exp_test
```

### 2. Using Pre-computed Coefficients
You can bypass the coefficients generation by providing the coefficient file path using the -f flag:
```bash
# Example: Running GeLU with pre-computed coefficients
./bin/examples/pke/gelu_test -f ../coeffs/gelucoeff.txt
```

## Plaintext Prototype
We provide a Python-based plaintext simulation tool in the project root directory: `fourier_extension.py`, which is used to prototype the Fourier Extension algorithm and find efficient parameters before running the homomorphic version.

### Examples & Usages
You can modify the configuration in the `if __name__ == "__main__":` block of `fourier_extension.py`:

- Function: Update `fstr` (e.g., `"x"`, `"exp(x)"`).
- Interval: Adjust `left` and `right`.
- Degree: Change `N` to the desired Fourier series order.

Run the script directly:
```bash
python fourier_extension.py
```
The script will output the recommended smoothness parameter, the estimated bit precision, and the generated coefficients:
```bash
ke: 33, Precision: 50.5166 bits
Fourier Coefficients (a_n, b_n): [(-3.885780586188048e-18+0j), (-1.8041124150158794e-18-0.625001603545843j), ....]
```
