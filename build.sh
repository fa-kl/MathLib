#!/usr/bin/env bash
set -e

# Colors for clarity
GREEN="\033[1;32m"
BLUE="\033[1;34m"
YELLOW="\033[1;33m"
RESET="\033[0m"

BUILD_DIR="build"
BUILD_TYPE="Release"
LIB_NAME="MathLib"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --coverage)
            ENABLE_COVERAGE="-DENABLE_COVERAGE=ON"
            shift
            ;;
        --install)
            INSTALL=true
            shift
            ;;
        --package)
            PACKAGE=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --debug      Build in Debug mode (default: Release)"
            echo "  --coverage   Enable code coverage"
            echo "  --install    Install libraries to install/ directory"
            echo "  --package    Create distribution packages"
            echo "  -h, --help   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Build Libraries (both static and shared)
# -----------------------------------------------------------------------------
echo -e "${GREEN}> Building MathLib (${BUILD_TYPE})...${RESET}"
echo -e "${BLUE}  • Static library: lib${LIB_NAME}.a${RESET}"
echo -e "${BLUE}  • Shared library: lib${LIB_NAME}.so${RESET}"

mkdir -p $BUILD_DIR && cd $BUILD_DIR

cmake .. \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_INSTALL_PREFIX=../install \
    ${ENABLE_COVERAGE}

make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# -----------------------------------------------------------------------------
# Show build results
# -----------------------------------------------------------------------------
echo -e "\n${GREEN}> Build completed! Generated libraries:${RESET}"
if [ -d "lib" ]; then
    ls -la lib/
else
    echo -e "${YELLOW}Warning: lib/ directory not found${RESET}"
fi

# -----------------------------------------------------------------------------
# Optional: Install
# -----------------------------------------------------------------------------
if [ "$INSTALL" = true ]; then
    echo -e "\n${GREEN}> Installing libraries...${RESET}"
    make install
    cd ..
    echo -e "${BLUE}Libraries installed to: install/${RESET}"
    echo -e "${BLUE}  • Headers: install/include/MathLib/${RESET}"
    echo -e "${BLUE}  • Libraries: install/lib/${RESET}"
    echo -e "${BLUE}  • CMake configs: install/lib/cmake/MathLib/${RESET}"
else
    cd ..
fi

# -----------------------------------------------------------------------------
# Optional: Package
# -----------------------------------------------------------------------------
if [ "$PACKAGE" = true ]; then
    echo -e "\n${GREEN}> Creating distribution packages...${RESET}"
    cd $BUILD_DIR
    cpack
    echo -e "${BLUE}Packages created in build/ directory:${RESET}"
    ls -la *.tar.gz *.zip *.deb *.rpm 2>/dev/null || echo "No packages found"
    cd ..
fi

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
echo -e "\n${GREEN}> All builds completed successfully!${RESET}"
echo -e "${YELLOW}Usage examples:${RESET}"
echo -e "  • Run tests: ./test.sh"
echo -e "  • Install: $0 --install"
echo -e "  • Create packages: $0 --package"
