#!/bin/bash

# ---------------------------------------------------------------------------------------
# coverage.sh - Generate code coverage report
#
# Usage:
#   ./coverage.sh          # Generate coverage report
#   ./coverage.sh clean    # Clean coverage data
# ---------------------------------------------------------------------------------------

set -e

# Check if lcov is installed
if ! command -v lcov &> /dev/null; then
    echo "Error: lcov is not installed."
    echo "Install it with: sudo apt-get install lcov"
    exit 1
fi

# Clean build directory if requested
if [ "$1" == "clean" ]; then
    echo "Cleaning coverage data..."
    rm -rf build/coverage_html
    rm -f build/coverage.info
    find build -name "*.gcda" -delete 2>/dev/null || true
    find build -name "*.gcno" -delete 2>/dev/null || true
    echo "Coverage data cleaned."
    exit 0
fi

# Create build directory if it doesn't exist
mkdir -p build
cd build

# Force clean reconfiguration for coverage mode
echo "Reconfiguring for coverage mode..."
rm -rf CMakeCache.txt CMakeFiles/

# Configure with coverage enabled
echo "Configuring CMake with coverage enabled..."
cmake -DENABLE_COVERAGE=ON ..

# Build the project
echo "Building project..."
cmake --build . -j$(nproc)

# Run coverage target
echo "Generating coverage report..."
cmake --build . --target coverage

# Display summary
echo ""
echo "=========================================="
echo "Coverage report generated successfully!"
echo "=========================================="
echo ""
echo "See:"
echo "  build/coverage_html/index.html"
echo ""

# Try to display a quick summary
if [ -f coverage.info ]; then
    echo "Coverage Summary:"
    lcov --summary coverage.info 2>&1 | grep -E "lines\.\.\.\.\.\.|functions\.\.\.\."
fi
