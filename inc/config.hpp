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

#include <stdint.h>

#include <array>
#include <cmath>

/// @brief Adds printing and string functionality to each type in MathLib.
#define __PRINT__

namespace mathlib
{

/// @brief A type for real values
using real32_t = float;
/// @brief A type for real values
using real64_t = double;
/// @brief A type for real values
using real_t = real64_t;
/// @brief A type for sizes and lengths
using size_t = std::size_t;
/// @brief A type for indices
using index_t = int64_t;

/// @brief A type for dimensions
typedef struct {
  size_t rows;
  size_t cols;
} dim_t;

/// @brief The value of π (pi)
const real_t PI = static_cast<real_t>(M_PI);
/// @brief The epsilon used for floating point comparisons
const real_t EPSILON = 1e-12;
/// @brief Infinity
#define INF INFINITY;

/// @brief Type trait to check if a type is a scalar numeric type
/// @details This includes all fundamental numeric types (int, float, double, etc.)
/// as well as the custom type aliases defined in this library (real_t, real32_t, real64_t)
template <typename T>
struct is_scalar_type {
  static constexpr bool value = std::is_arithmetic<T>::value;
};

/// @brief Helper variable template for is_scalar_type
template <typename T>
constexpr bool is_scalar_type_v = is_scalar_type<T>::value;

}  // namespace mathlib