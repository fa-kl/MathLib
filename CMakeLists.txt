# ---------------------------------------------------------------------------------------
#
# @file: CMakeLists.txt
#
# @brief: The library's main CMakeLists.txt file
#
# @author: fakl
# @date: October 2025
#
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# Setup & Settings
# ---------------------------------------------------------------------------------------
set(LIB_NAME "MathLib")
set(LIB_VERSION "1.0.0")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


cmake_minimum_required(VERSION 3.20)
project(${LIB_NAME} VERSION ${LIB_VERSION} LANGUAGES CXX)

# ---------------------------------------------------------------------------------------
# Coverage Options
# ---------------------------------------------------------------------------------------
option(ENABLE_COVERAGE "Enable code coverage" OFF)

# ---------------------------------------------------------------------------------------
# Warning Settings
# ---------------------------------------------------------------------------------------
# Enabling all possible warnings
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # GCC/Clang: Enable extra warnings
    add_compile_options(
        -Wall                   # Enable common warnings
        -Wextra                 # Enable extra warnings
        -Wpedantic              # Enable strict compliance warnings
        -Werror                 # Treat warnings as errors
        -Wshadow                # Warn about variable shadowing
        -Wformat=2              # Check for format string issues (use `-Wformat=2` for stricter)
        -Wconversion            # Warn on implicit type conversions
        -Wsign-conversion       # Warn about signed/unsigned comparisons
        -Wuninitialized         # Warn if a variable is used before initialization
        -Wno-unused-parameter   # Optional: disable unused parameter warnings if needed
    )
elseif(MSVC)
    # MSVC: Enable as many warnings as possible
    add_compile_options(
        /W4               # Enable high warning level
        /WX               # Treat warnings as errors
        /wd4800           # Ignore warning about 'forcing value to bool' (optional)
    )
elseif (APPLE)
    # Apple-specific flags (Clang)
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Werror
        -Wshadow
        -Wconversion
        -Wsign-conversion
        -Wuninitialized
    )
endif()

# ---------------------------------------------------------------------------------------
# Coverage Flags
# ---------------------------------------------------------------------------------------
if(ENABLE_COVERAGE AND (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    message(STATUS "Code coverage enabled")
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
endif()

# ---------------------------------------------------------------------------------------
# Library Sources & Headers
# ---------------------------------------------------------------------------------------
file(GLOB LIB_SOURCES src/*.cpp)

# Create both static and shared libraries
add_library("${LIB_NAME}_static" STATIC ${LIB_SOURCES})
add_library("${LIB_NAME}_shared" SHARED ${LIB_SOURCES})

# Set output names to avoid conflicts
set_target_properties("${LIB_NAME}_static" PROPERTIES
    OUTPUT_NAME "${LIB_NAME}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

set_target_properties("${LIB_NAME}_shared" PROPERTIES
    OUTPUT_NAME "${LIB_NAME}"
    VERSION ${LIB_VERSION}
    SOVERSION 1
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Include directories for both libraries
target_include_directories("${LIB_NAME}_static" PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)

target_include_directories("${LIB_NAME}_shared" PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)

# Create an alias for backward compatibility (defaults to static)
add_library("${LIB_NAME}" ALIAS "${LIB_NAME}_static")

# ---------------------------------------------------------------------------------------
# Main Executable for Debugging
# ---------------------------------------------------------------------------------------
set(MAIN_TARGET "${LIB_NAME}Main")
add_executable(${MAIN_TARGET} main.cpp)
target_link_libraries(${MAIN_TARGET} PRIVATE "${LIB_NAME}")

# ---------------------------------------------------------------------------------------
# Testing
# ---------------------------------------------------------------------------------------
include(CTest)
enable_testing()
include(FetchContent)

# Configure GoogleTest to not install itself
set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG origin/main)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

set(TEST_TARGET "${LIB_NAME}Tests")
file(GLOB TEST_SOURCES test/*.cpp)
add_executable(${TEST_TARGET} ${TEST_SOURCES})
target_link_libraries(${TEST_TARGET} PRIVATE "${LIB_NAME}" gtest_main)

# Disable self-move warning for tests (we intentionally test self-move behavior)
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${TEST_TARGET} PRIVATE -Wno-self-move)
elseif(MSVC)
    # MSVC doesn't have this specific warning
endif()

gtest_discover_tests(${TEST_TARGET})

# ---------------------------------------------------------------------------------------
# Coverage Target
# ---------------------------------------------------------------------------------------
if(ENABLE_COVERAGE)
    find_program(LCOV lcov)
    find_program(GENHTML genhtml)
    
    if(LCOV AND GENHTML)
        add_custom_target(coverage
            COMMAND ${CMAKE_COMMAND} -E echo "Running tests..."
            COMMAND ${TEST_TARGET} || true
            
            COMMAND ${CMAKE_COMMAND} -E echo "Capturing coverage data..."
            COMMAND ${LCOV} --capture --directory . --output-file coverage.info --ignore-errors mismatch,gcov --rc geninfo_unexecuted_blocks=1
            
            COMMAND ${CMAKE_COMMAND} -E echo "Filtering coverage data..."
            COMMAND ${LCOV} --remove coverage.info 
                '/usr/*' 
                '*/test/*' 
                '*/googletest-*/*'
                '*/build/_deps/*'
                --output-file coverage.info --ignore-errors unused
            
            COMMAND ${CMAKE_COMMAND} -E echo "Generating HTML report..."
            COMMAND ${GENHTML} coverage.info --output-directory coverage_html --ignore-errors source
            
            COMMAND ${CMAKE_COMMAND} -E echo "Coverage report generated in coverage_html/index.html"
            
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating code coverage report"
        )
        
        add_custom_target(coverage-clean
            COMMAND ${CMAKE_COMMAND} -E remove_directory coverage_html
            COMMAND ${CMAKE_COMMAND} -E remove coverage.info
            COMMAND find . -name "*.gcda" -delete
            COMMAND find . -name "*.gcno" -delete
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Cleaning coverage data"
        )
    else()
        message(WARNING "lcov and/or genhtml not found. Coverage target not available.")
    endif()
endif()

# ---------------------------------------------------------------------------------------
# Installation Configuration
# ---------------------------------------------------------------------------------------
include(GNUInstallDirs)

# Install libraries
install(TARGETS "${LIB_NAME}_static" "${LIB_NAME}_shared"
    EXPORT "${LIB_NAME}Targets"
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install headers
install(DIRECTORY inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${LIB_NAME}
    FILES_MATCHING PATTERN "*.hpp"
)

# Install CMake config files
install(EXPORT "${LIB_NAME}Targets"
    FILE "${LIB_NAME}Targets.cmake"
    NAMESPACE "${LIB_NAME}::"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
    VERSION ${LIB_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${LIB_NAME}Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
)

# ---------------------------------------------------------------------------------------
# CPack Configuration for Release Packaging
# ---------------------------------------------------------------------------------------
include(CPack)

set(CPACK_PACKAGE_NAME "${LIB_NAME}")
set(CPACK_PACKAGE_VERSION "${LIB_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MathLib - A comprehensive C++ mathematics library")
set(CPACK_PACKAGE_VENDOR "fakl")
set(CPACK_PACKAGE_CONTACT "your-email@example.com")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

# Package formats
set(CPACK_GENERATOR "TGZ;ZIP")
if(UNIX AND NOT APPLE)
    list(APPEND CPACK_GENERATOR "DEB;RPM")
endif()

# DEB package configuration
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "fakl")
set(CPACK_DEBIAN_PACKAGE_SECTION "libs")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

# RPM package configuration
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
set(CPACK_RPM_PACKAGE_LICENSE "See LICENSE file")

# Include both static and shared libraries in packages
set(CPACK_COMPONENTS_ALL libraries headers cmake)
