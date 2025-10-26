/*****************************************************************************************
 * @file: compare.hpp
 *
 * @brief: Comparison functions
 *
 * @details: This file provides various floating point comparison functions.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include "config.hpp"

namespace mathlib {

/// @brief Compare to real (floating point) values for equality
/// @param value1 First real value
/// @param value2 Second real value
/// @param epsilon Comparison tolerance
/// @returns True if both values are equal within the given tolerance
bool isFuzzyEqual(real_t value1, real_t value2, real_t epsilon = EPSILON);

/// @brief Check whether the value1 is greater or equal than value2
/// @param value1 First real value
/// @param value2 Second real value
/// @param epsilon Comparison tolerance
/// @returns True if the first value is greater or equal than the second value
/// within the given tolerance
bool isFuzzyGreater(real_t value1, real_t value2, real_t epsilon = EPSILON);

/// @brief Check whether the value1 is strictly smaller than value2
/// @param value1 First real value
/// @param value2 Second real value
/// @param epsilon Comparison tolerance
/// @returns True if the first value is smaller or equal than the second value
/// within the given tolerance
bool isFuzzySmaller(real_t value1, real_t value2, real_t epsilon = EPSILON);

/// @brief Check whether the value1 is strictly greater than value2
/// @param value1 First real value
/// @param value2 Second real value
/// @param epsilon Comparison tolerance
/// @returns True if the first value is strictly greater than the second value
/// within the given tolerance
bool isStrictFuzzyGreater(real_t value1, real_t value2,
                          real_t epsilon = EPSILON);

/// @brief Check whether the value1 is strictly smaller than value2
/// @param value1 First real value
/// @param value2 Second real value
/// @param epsilon Comparison tolerance
/// @returns True if the first value is strictly smaller than the second value
/// within the given tolerance
bool isStrictFuzzySmaller(real_t value1, real_t value2,
                          real_t epsilon = EPSILON);

} // namespace mathlib
