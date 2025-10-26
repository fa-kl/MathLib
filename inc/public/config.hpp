/*****************************************************************************************
 * @file: config.hpp
 *
 * @brief: MathLib's configuration file
 *
 * @details: This configures basic types and other settings.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <cmath>
#include <stdint.h>

/// @brief Adds printing and string functionality to each type in MathLib.
#define __PRINT__

namespace mathlib {

/// @brief A type for real values
using real_t = double;
/// @brief A type for sizes and lengths
using size_t = uint32_t;
/// @brief A type for indices
using index_t = int64_t;
/// @brief The value of π (pi)
const real_t PI = M_PI;
/// @brief The epsilon used for floating point comparisons
const real_t EPSILON = 1e-9;

} // namespace mathlib