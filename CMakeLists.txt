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
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


cmake_minimum_required(VERSION 3.20)
project(${LIB_NAME})

# ---------------------------------------------------------------------------------------
# Library Sources & Headers
# ---------------------------------------------------------------------------------------
# Static Library
file(GLOB LIB_SOURCES src/*.cpp)
add_library("${LIB_NAME}" STATIC ${LIB_SOURCES})
target_include_directories("${LIB_NAME}" PUBLIC inc/public PRIVATE inc/private)

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
gtest_discover_tests(${TEST_TARGET})