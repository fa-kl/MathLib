/*****************************************************************************************
 * @file: rounding.hpp
 *
 * @brief: Rounding capabilities.
 *
 * @details: This file declares rounding methods.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <math.h>

#include "MathLibError.hpp"
#include "config.hpp"

namespace mathlib
{

typedef enum { NEAREST, CEIL, FLOOR } RoundingMethod_t;

//------------------------------------------------------------------------------
// Integer rounding (no decimal precision, but with method selection)
//------------------------------------------------------------------------------
inline int64_t round(const real32_t& value, RoundingMethod_t method = NEAREST)
{
  switch (method) {
    case NEAREST:
      return static_cast<int64_t>(std::round(value));
    case CEIL:
      return static_cast<int64_t>(std::ceil(value));
    case FLOOR:
      return static_cast<int64_t>(std::floor(value));
    default:
      throw MathLibError("Invalid rounding method");
  }
}

//------------------------------------------------------------------------------
// Floating-point rounding (float)
//------------------------------------------------------------------------------
inline real32_t round(const real32_t& value, const size_t& decimal, RoundingMethod_t method = NEAREST)
{
  real32_t scale = static_cast<real32_t>(std::pow(10.0f, static_cast<real32_t>(decimal)));
  switch (method) {
    case NEAREST:
      return std::roundf(value * scale) / scale;
    case CEIL:
      return std::ceil(value * scale) / scale;
    case FLOOR:
      return std::floor(value * scale) / scale;
    default:
      throw MathLibError("Invalid rounding method");
  }
}

//------------------------------------------------------------------------------
// Floating-point rounding (double)
//------------------------------------------------------------------------------
inline real64_t round(const real64_t& value, const size_t& decimal, RoundingMethod_t method = NEAREST)
{
  real64_t scale = static_cast<real64_t>(std::pow(10.0, static_cast<real64_t>(decimal)));
  switch (method) {
    case NEAREST:
      return std::round(value * scale) / scale;
    case CEIL:
      return std::ceil(value * scale) / scale;
    case FLOOR:
      return std::floor(value * scale) / scale;
    default:
      throw MathLibError("Invalid rounding method");
  }
}

//------------------------------------------------------------------------------
// Convenience wrappers (float)
//------------------------------------------------------------------------------
inline real32_t nearest(const real32_t& value, const size_t& decimal)
{
  return round(value, decimal, NEAREST);
}

inline real32_t ceil(const real32_t& value, const size_t& decimal)
{
  return round(value, decimal, CEIL);
}

inline real32_t floor(const real32_t& value, const size_t& decimal)
{
  return round(value, decimal, FLOOR);
}

//------------------------------------------------------------------------------
// Convenience wrappers (double)
//------------------------------------------------------------------------------
inline real64_t nearest(const real64_t& value, const size_t& decimal)
{
  return round(value, decimal, NEAREST);
}

inline real64_t ceil(const real64_t& value, const size_t& decimal)
{
  return round(value, decimal, CEIL);
}

inline real64_t floor(const real64_t& value, const size_t& decimal)
{
  return round(value, decimal, FLOOR);
}

}  // namespace mathlib