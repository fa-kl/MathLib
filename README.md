# MathLib
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/fa-kl/MathLib/releases)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![C++](https://img.shields.io/badge/C++-23-blue.svg)](https://en.cppreference.com/w/cpp/23)

A modern C++23 library for mathematical operations including complex numbers, vectors, matrices, and linear algebra. 

## Features
  - Floating-point comparison
  - Unit conversions
  - Complex numbers
  - n-dimensional vectors
  - Matrices of arbitrary size
    
### Floating-Point Comparison
  - Comparison for equality: ```isFuzzyEqual()```
  - Check for strictly greater: ```isStrictFuzzyGreater()```
  - Check for greater or equal: ```isFuzzyGreater()```
  - Check for strictly smaller: ```isStrictFuzzySmaller()```
  - Check for smaller or equal: ```isFuzzySmaller()```
    
### Unit Conversions
  - Radiants to degrees: ```rad2deg()```
  - Degrees to radiants: ```deg2rad()```
    
### Complex Numbers (```class Complex```)
  - Getters: ```real()``` and ```imag()```
  - Arithmetic operations: ```+```, ```-```, ```*```, ```/```, ```+=```, ```-=```, ```*=```, ```/=```
    - ```Complex``` & ```Complex```
    - ```Complex``` & ```real_t```
    - ```real_t``` & ```Complex```
  - Magnitude: ```abs()```
  - Squared magnitude: ```abs2()```
  - Argument: ```arg()```
  - Conjungate: ```conj()```
  - Exponential map: ```exp()```
  - Logarithmic map: ```log()```
  - Trigonometric maps: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Power: ```pow()```
  - Square root: ```sqrt()```
  - Conversion to ```std::string```: ```to_string()```
  - Serialization and printing to the default output stream: ```<<``` and ```print()```

### n-Dimensional Vectors (```Vector```)
  - Constructors:
    - Empty vector: ```Vector()```
    - n-dimensional vector: ```Vector(n)```
    - Copy & move constructors and operators
    - Vector of all ones: ```ones(n)```
    - Vector of all zeros: ```zeros(n)```
  - Element access
    - 0-based: ```operator[]```
    - 1-based: ```operator()``` allows positive & negative indices
  - Information acess
    - Number of elements: ```length()```
    - Dimension: ```size()```
    - Number of rows: ```rows()```
    - Number of columns: ```cols()```
  - Element-wise addition and subtraction
    - ```Vector``` & ```Vector```: ```+```, ```-```
    - ```Vector``` & scalar: ```+```, ```-```
    - scalar & ```Vector```: ```+```, ```-```
  - Element-wise multiplicaiton and division
    - Vector-Vector: ```elmul()```, ```eldiv()```
    - Vector-Scalar: ```*```, ```/```
    - Scalar-Vector: ```*```, ```/```
  - Dot product: ```*```, ```dot()```
  - Cross product: ```cross()```
  - Outer product: ```outer()```
  - Transpose: ```transpose()```
  - Norm: ```norm()```
  - Normalization: ```normalize()```
  - Sum of elements: ```sum()```
  - Difference between vector elements: ```diff()```
  - Element-wise rounding: ```round()```
  - Element-wise absolute values: ```abs()```
  - Maximum and minimum: ```max()```, ```min()```
  - Mean, standard deviation, variance: ```mean()```, ```std()```, ```var()```
  - Exponential and logarithm: ```exp()```, ```log()```
  - Element-wise square root and pow: ```sqrt()```, ```pow()```
  - Element-wise trigonometric functions: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Element-wise logical operators: ```&&```, ```||```, ```^```
  - Element-wise comparison via operators: ```==```, ```!=```, ```>```, ```>=```, ```<```, ```<=```
  - Element-wise floating-point comparison via methods: ```isFuzzyEqual()```, ```isStrictFuzzyGreater()```, ```isFuzzyGreater()```, ```isStrictFuzzySmaller()```, ```isFuzzySmaller()```
  - Conversion to ```std::string``` and printing via ```<iostream>```
  - and more ...
    
### Matrices (```Matrix```)  
  - Constructors:
    - Empty matrix: ```Matrix()```
    - (m, n)-dimensional matrix: ```Matrix(m, n)```
    - Copy & move constructors and operators
    - Identity matrix: ```eye(n)```
    - Matrix of all ones: ```ones(m, n)```
    - Matrix of all zeros: ```zeros(m, n)```
    - Diagonal matrix: ```diag(Vector)``` or ```diag({a, b, ...})```
  - Element access
    - 0-based: ```operator[]``` (linear access only)
    - 1-based: ```operator()``` (positive & negative linear and positive & negative 2D indices)
  - Information acess
    - Number of elements: ```length()```
    - Dimension: ```size()```
    - Number of rows: ```rows()```
    - Number of columns: ```cols()```
    - Rank: ```rank()```
  - Element-wise addition and subtraction
    - ```Matrix``` & ```Matrix```: ```+```, ```-```
    - ```Matrix``` & scalar: ```+```, ```-```
    - scalar & ```Matrix```: ```+```, ```-```
  - Element-wise multiplicaiton
    - Matrix-Matrix: ```*```
    - Matrix-Vector: ```*```
    - Vector-Matrix: ```*```
    - Matrix-Scalar: ```*```
    - Scalar-Matrix: ```*```
  - Element-wise division
    - Matrix-scalar: ```/```
  - Dot product: ```*```, ```dot()```, ```colDot()```, ```rowDot()```
  - Outer product: ```outer()```
  - Transpose: ```transpose()```
  - Trace:```trace()``` 
  - Matrix norms: ```norm()```, ```vecnorm()```
  - Orthogonalization: ```orth()```
  - QR decomposition: ```qr()```
  - Eigenvalues and -vectors: ```eig()```
  - SVD: ```svd()```
  - Normalization: ```normalize()``` (rows or columns)
  - Sum of elements: ```sum()``` (all, rows-wise, column-wise)
  - Element-wise rounding: ```round()```
  - Element-wise absolute values: ```abs()```
  - Maximum and minimum: ```max()```, ```min()``` (all, rows-wise, column-wise)
  - Mean, standard deviation, variance: ```mean()```, ```std()```, ```var()```  (all, rows-wise, column-wise)
  - Exponential and logarithm: ```exp()```, ```log()```
  - Element-wise square root and pow: ```sqrt()```, ```elpow()```
  - Element-wise trigonometric functions: ```sin()```, ```cos()```, ```tan()```, ```asin()```, ```acos()```, ```atan()```
  - Element-wise logical operators: ```&&```, ```||```, ```^```
  - Element-wise comparison via operators: ```==```, ```!=```, ```>```, ```>=```, ```<```, ```<=```
  - Element-wise floating-point comparison via methods: ```isFuzzyEqual()```, ```isStrictFuzzyGreater()```, ```isFuzzyGreater()```, ```isStrictFuzzySmaller()```, ```isFuzzySmaller()```
  - Conversion to ```std::string``` and printing via ```<iostream>```
  - and more ...

## Requirements
- **C++23** compatible compiler:
  - GCC 12 or higher
  - Clang 16 or higher
  - MSVC 2022 or higher (Visual Studio 2022)
- **CMake 3.20** or higher

## Installation

### Option 1: Pre-built Libraries (Recommended)

Download pre-compiled libraries from the [GitHub Releases](https://github.com/fa-kl/MathLib/releases) page.

#### 1. Download the appropriate package
- **Linux**: Download `mathlib-linux.tar.gz`
- **Windows**: Download `mathlib-windows.zip`
- **macOS**: Download `mathlib-macos.tar.gz`

#### 2. Extract the package
```bash
# Linux/macOS
tar -xzf mathlib-linux.tar.gz
cd mathlib-linux/

# Windows (PowerShell)
Expand-Archive mathlib-windows.zip
cd mathlib-windows/
```

#### 3. Choose your library type

Each package contains both static and shared libraries:
- **Static library**: `libMathLib.a` (Unix) or `MathLib.lib` (Windows)
- **Shared library**: `libMathLib.so` (Linux), `libMathLib.dylib` (macOS), or `MathLib.dll` (Windows)

#### 4. Integrate with your project

##### Using CMake (Recommended)
```cmake
# Method 1: Use the provided CMake config
find_package(MathLib REQUIRED PATHS /path/to/extracted/mathlib/cmake)

# Link static library (default)
target_link_libraries(your_target PRIVATE MathLib::MathLib_static)

# Or link shared library
target_link_libraries(your_target PRIVATE MathLib::MathLib_shared)
```

##### Manual Integration
```cmake
# Set paths
set(MATHLIB_ROOT "/path/to/extracted/mathlib")
set(MATHLIB_INCLUDE_DIR "${MATHLIB_ROOT}/include/MathLib")

# Include headers
target_include_directories(your_target PRIVATE ${MATHLIB_INCLUDE_DIR})

# Link static library
target_link_libraries(your_target PRIVATE "${MATHLIB_ROOT}/lib/libMathLib.a")

# Or link shared library (Linux)
target_link_libraries(your_target PRIVATE "${MATHLIB_ROOT}/lib/libMathLib.so")
```

##### Using pkg-config (Linux/macOS)
```bash
# Copy .pc file to pkg-config directory or set PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/path/to/extracted/mathlib/lib/pkgconfig:$PKG_CONFIG_PATH

# In your Makefile or build script
CFLAGS += $(shell pkg-config --cflags mathlib)
LDFLAGS += $(shell pkg-config --libs mathlib)
```

### Option 2: Build from Source

#### 1. Clone the repository
```bash
git clone https://github.com/fa-kl/MathLib.git
cd MathLib
```

#### 2. Build both static and shared libraries
```bash
./build.sh --install --package
```

This creates:
- Libraries in `build/lib/`
- Headers and CMake configs in `install/`
- Distribution packages in `build/`

#### 3. Manual build with CMake
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc)
make install
```

#### 4. Link to your project
```cmake
# Add MathLib as subdirectory
add_subdirectory(path/to/MathLib)

# Link static library (default)
target_link_libraries(your_target PRIVATE MathLib::MathLib_static)

# Or link shared library
target_link_libraries(your_target PRIVATE MathLib::MathLib_shared)
```

### Option 3: Git Submodule

Add MathLib as a submodule to your project:
```bash
git submodule add https://github.com/fa-kl/MathLib.git external/MathLib
git submodule update --init --recursive
```

Then in your `CMakeLists.txt`:
```cmake
add_subdirectory(external/MathLib)
# Use static library by default, or specify explicitly
target_link_libraries(your_target PRIVATE MathLib::MathLib_static)
# target_link_libraries(your_target PRIVATE MathLib::MathLib_shared)
```

### Option 4: Header-only Integration

For simple projects, copy the `inc/` and `src/` directories:

```cpp
#include "mathlib.hpp"  // Main header that includes all components
```

### Library Type Comparison

| Feature | Static Library | Shared Library |
|---------|---------------|----------------|
| **File size** | Larger executable | Smaller executable |
| **Dependencies** | Self-contained | Requires .so/.dll/.dylib |
| **Loading time** | Faster startup | Slightly slower startup |
| **Memory usage** | Higher (if multiple programs use it) | Lower (shared between programs) |
| **Distribution** | Easier (single file) | Requires library file |
| **Updates** | Need to recompile | Can update library separately |

**Recommendation**: Use static libraries for simpler distribution, shared libraries for system-wide installation or when multiple applications use MathLib.

## Quick Start Example

```cpp
#include "mathlib.hpp"
using namespace mathlib;

int main() {
    // Complex numbers
    Complex z1(3.0, 4.0);
    Complex z2 = exp(z1);
    
    // Vectors
    Vector<double> v1 = {1.0, 2.0, 3.0};
    Vector<double> v2 = {4.0, 5.0, 6.0};
    double dot = v1 * v2;  // Dot product
    
    // Matrices
    Matrix<double> A = {{1.0, 2.0}, 
                        {3.0, 4.0}};
    Matrix<double> B = A.transpose();
    auto qr = mathlib::qr(A);  // QR decomposition
    
    return 0;
}
```

## Running Tests

```bash
cd build
ctest
# or
./MathLibTests
```

## Running with Coverage

```bash
mkdir build
cd build
cmake -DENABLE_COVERAGE=ON ..
make
make coverage
# Open coverage_html/index.html in browser
```

Or use the provided script:
```bash
./coverage.sh
```

## Documentation

Full API documentation is available at: **[https://fa-kl.github.io/MathLib](https://fa-kl.github.io/MathLib)**

The documentation is generated from the source code using Doxygen. To generate it locally:
```bash
doxygen Doxyfile
```

## Project Structure

```
MathLib/
├── inc/            # Header files
│   ├── mathlib.hpp # Main include file
│   ├── Complex.hpp
│   ├── Vector.hpp
│   ├── Matrix.hpp
│   └── ...
├── src/            # Implementation files
├── test/           # Unit tests (GoogleTest)
├── build/          # Build directory (generated)
└── CMakeLists.txt  # CMake configuration
```

## Contributing

This is primarily a learning project to improve C++ skills. If you want to improve your C++ skills or have any tips, feel free to:
- Open an issue
- Submit a pull request
- Contact me for discussions

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and changes.
