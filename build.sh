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
        --release)
            RELEASE=true
            RELEASE_VERSION="$2"
            if [[ -z "$RELEASE_VERSION" ]]; then
                echo "Error: --release requires a version (e.g., --release v1.0.2)"
                exit 1
            fi
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --debug      Build in Debug mode (default: Release)"
            echo "  --coverage   Enable code coverage"
            echo "  --install    Install libraries to install/ directory"
            echo "  --package    Create distribution packages"
            echo "  --release <version>  Create release package (e.g., --release v1.0.2)"
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
# Optional: Create Release Package
# -----------------------------------------------------------------------------
if [ "$RELEASE" = true ]; then
    # Ensure we have built and installed first
    if [ ! -d "install" ]; then
        echo -e "\n${YELLOW}> Release requires installation, building and installing first...${RESET}"
        INSTALL=true
        cd $BUILD_DIR
        make install
        cd ..
    fi
    
    echo -e "\n${GREEN}> Creating release package: ${RELEASE_VERSION}${RESET}"
    
    # Create release directory structure
    RELEASE_DIR="release/${RELEASE_VERSION}"
    mkdir -p "${RELEASE_DIR}/lib" "${RELEASE_DIR}/include/MathLib" "${RELEASE_DIR}/cmake"
    
    # Copy libraries
    cp install/lib/libMathLib.a install/lib/libMathLib.so* "${RELEASE_DIR}/lib/"
    
    # Copy headers
    cp -r install/include/MathLib/* "${RELEASE_DIR}/include/MathLib/"
    
    # Copy CMake configuration
    cp -r install/lib/cmake/MathLib "${RELEASE_DIR}/cmake/"
    
    # Copy documentation
    cp README.md LICENSE "${RELEASE_DIR}/"
    
    # Copy installation guide
    cp release/INSTALL.md "${RELEASE_DIR}/"
    
    # Create compressed archive
    ARCHIVE_NAME="release/mathlib-${RELEASE_VERSION}-linux-x64.tar.gz"
    tar -czf "${ARCHIVE_NAME}" -C "${RELEASE_DIR}" .
    
    echo -e "${BLUE}Release package created:${RESET}"
    echo -e "${BLUE}  • Directory: ${RELEASE_DIR}${RESET}"
    echo -e "${BLUE}  • Archive: ${ARCHIVE_NAME} ($(du -h ${ARCHIVE_NAME} | cut -f1))${RESET}"
    echo -e "${BLUE}  • Contents: libraries, headers, CMake configs, docs${RESET}"
fi

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
echo -e "\n${GREEN}> All builds completed successfully!${RESET}"
echo -e "${YELLOW}Usage examples:${RESET}"
echo -e "  • Run tests: ./test.sh"
echo -e "  • Install: $0 --install"
echo -e "  • Create packages: $0 --package"
echo -e "  • Create release: $0 --release v1.0.2"
