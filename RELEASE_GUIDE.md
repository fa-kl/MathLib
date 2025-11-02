# Creating the v1.0.0 Release

Follow these steps to create your first release on GitHub.

## Pre-Release Checklist ✅

- [x] All tests passing
- [x] Version number updated in CMakeLists.txt (1.0.0)
- [x] CHANGELOG.md created with release notes
- [x] README.md updated with installation instructions
- [x] Documentation setup (Doxygen + GitHub Pages)
- [x] .gitignore updated
- [x] License file exists (MIT)

## Step-by-Step Release Process

### Step 1: Commit and Push All Changes

```bash
# Check status
git status

# Add all new files
git add .

# Commit
git commit -m "Release v1.0.0: Initial stable release with full documentation"

# Push to main
git push origin main
```

### Step 2: Create a Git Tag

```bash
# Create an annotated tag
git tag -a v1.0.0 -m "MathLib v1.0.0 - Initial Release"

# Push the tag to GitHub
git push origin v1.0.0
```

### Step 3: Create Release on GitHub

1. Go to https://github.com/fa-kl/MathLib
2. Click **Releases** (right sidebar, under "About")
3. Click **"Create a new release"** or **"Draft a new release"**
4. Fill in the release form:

#### Choose a tag
   - Select: **v1.0.0** (the tag you just pushed)
   - Or type `v1.0.0` to create on publish

#### Release title
```
MathLib v1.0.0 - Initial Release
```

#### Release description
```markdown
# 🎉 MathLib v1.0.0 - Initial Release

First stable release of **MathLib**, a modern C++23 mathematical library for complex numbers, vectors, matrices, and linear algebra operations.

## ✨ Features

### Core Components
- **Complex Numbers**: Full arithmetic with trigonometric, exponential, and logarithmic functions
- **N-Dimensional Vectors**: Comprehensive vector operations including dot/cross products, norms, and statistical functions
- **Matrices**: Full matrix algebra with QR decomposition, eigenvalue computation, and rank analysis

### Mathematical Operations
- Element-wise arithmetic and logical operations
- Fuzzy floating-point comparisons for numerical stability
- Trigonometric functions (sin, cos, tan, asin, acos, atan)
- Exponential, logarithm, power, and square root operations
- Statistical functions (mean, std, var, max, min)
- Linear algebra (QR decomposition, eigenvalues/eigenvectors, orthogonalization)

### Development Features
- Comprehensive test suite using GoogleTest
- CMake build system with coverage support
- Complete Doxygen documentation
- Cross-platform compiler support (GCC, Clang, MSVC)

## 📋 Requirements

- **C++23** compatible compiler (GCC 12+, Clang 16+, MSVC 2022+)
- **CMake 3.20** or higher

## 📚 Documentation

Full API documentation: https://fa-kl.github.io/MathLib

## 🚀 Quick Start

```cpp
#include "mathlib.hpp"
using namespace mathlib;

int main() {
    // Complex numbers
    Complex z(3.0, 4.0);
    auto magnitude = abs(z);
    
    // Vectors
    Vector<double> v = {1.0, 2.0, 3.0};
    auto normalized = normalize(v);
    
    // Matrices
    Matrix<double> A = {{1.0, 2.0}, {3.0, 4.0}};
    auto qr_result = qr(A);
    
    return 0;
}
```

## 📦 Installation

See the [README](https://github.com/fa-kl/MathLib/blob/main/README.md) for detailed installation instructions.

**Quick CMake Integration:**
```cmake
add_subdirectory(path/to/MathLib)
target_link_libraries(your_target PRIVATE MathLib)
```

## 🔄 What's Next?

Future plans include:
- SVD (Singular Value Decomposition) implementation
- Additional matrix decomposition methods
- Performance optimizations
- More examples and tutorials

## 📄 License

This project is licensed under the MIT License.

---

**Full Changelog**: https://github.com/fa-kl/MathLib/blob/main/CHANGELOG.md
```

5. **Set as the latest release**: ✅ Check this box
6. Click **"Publish release"**

### Step 4: Verify the Release

1. Check that the release appears at: https://github.com/fa-kl/MathLib/releases
2. Verify the tag exists: https://github.com/fa-kl/MathLib/tags
3. Check that documentation deploys: https://fa-kl.github.io/MathLib (may take 5-10 minutes)

## Post-Release Actions

### Share Your Release! 🎊

You can now share your library:
- The release link: `https://github.com/fa-kl/MathLib/releases/tag/v1.0.0`
- The documentation: `https://fa-kl.github.io/MathLib`

### Add Release Badge to README (Optional)

The badges are already in your README:
- Version badge will show v1.0.0
- License badge shows MIT
- C++ version badge shows C++23

### Future Releases

For subsequent releases (v1.0.1, v1.1.0, etc.):
1. Update version in `CMakeLists.txt`
2. Update `CHANGELOG.md` with new changes
3. Update `Doxyfile` PROJECT_NUMBER
4. Commit, tag, and create release following the same process

## Semantic Versioning Guide

- **MAJOR** (v2.0.0): Breaking API changes
- **MINOR** (v1.1.0): New features, backward compatible
- **PATCH** (v1.0.1): Bug fixes, backward compatible

---

**Congratulations on your first release! 🎉**
