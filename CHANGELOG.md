# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-02

### Added
- Initial release of MathLib
- Complex number class with full arithmetic and mathematical operations
- N-dimensional vector class with extensive operations:
  - Arithmetic operations (addition, subtraction, multiplication, division)
  - Dot product, cross product, outer product
  - Norms and normalization
  - Statistical functions (mean, std, var, max, min)
  - Trigonometric functions (sin, cos, tan, asin, acos, atan)
  - Exponential, logarithm, power, square root
  - Element-wise operations
  - Fuzzy floating-point comparisons
- Matrix class with comprehensive linear algebra support:
  - Matrix arithmetic (addition, subtraction, multiplication)
  - Matrix-vector operations
  - Transpose, trace, determinant (via rank)
  - QR decomposition
  - Eigenvalue and eigenvector computation
  - Rank computation and orthogonalization
  - Statistical operations (mean, std, var) with dimension support
  - Matrix norms
  - Element-wise operations and comparisons
- Floating-point comparison utilities
- Unit conversion utilities (degrees/radians)
- Comprehensive test suite using GoogleTest
- CMake build system with:
  - Static library compilation
  - Test integration
  - Code coverage support
  - Cross-platform compiler warning configurations
- Complete Doxygen documentation comments
- MIT License

### Requirements
- C++23 compatible compiler
- CMake 3.20 or higher

### Known Limitations
- SVD implementation is mentioned in README but not yet implemented
- Documentation requires manual Doxygen generation

[1.0.0]: https://github.com/fa-kl/MathLib/releases/tag/v1.0.0
