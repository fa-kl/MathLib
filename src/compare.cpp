/*****************************************************************************************
 * @file: compare.cpp
 *
 * @brief: Comparison functions
 *
 * @details: This file provides various floating point comparison functions.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "compare.hpp"
#include <cmath>

namespace mathlib {

bool isFuzzyEqual(real_t value1, real_t value2, real_t epsilon) {
  return std::abs(value1 - value2) < epsilon;
}

bool isFuzzyGreater(real_t value1, real_t value2, real_t epsilon) {
  return value1 + epsilon > value2;
}

bool isFuzzySmaller(real_t value1, real_t value2, real_t epsilon) {
  return value1 < epsilon + value2;
}

bool isStrictFuzzyGreater(real_t value1, real_t value2, real_t epsilon) {
  return value1 > epsilon + value2;
}

bool isStrictFuzzySmaller(real_t value1, real_t value2, real_t epsilon) {
  return value1 + epsilon < value2;
}

} // namespace mathlib