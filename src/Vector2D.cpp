/*****************************************************************************************
 * @file: Vector2D.hpp
 *
 * @brief: A class for 2D vectors.
 *
 * @details: This file implements a class Vector2D for handling 2D vectors.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "Vector2D.hpp"
#include "DivByZeroError.hpp"
#include "IndexOutOfRangeError.hpp"
#include "MathLibError.hpp"
#include "compare.hpp"
#include "macros.h"

#include <cmath>
#include <exception>
#ifdef __PRINT__
#include <iostream>
#endif

namespace mathlib {

/// @brief Create a 2D vector of all zeros
Vector2D::Vector2D() : m_x(0.0), m_y(0.0) {}

/// @brief Create a 2D vector
/// @param x The vector's x component
/// @param y The vector's y component
Vector2D::Vector2D(const real_t x, const real_t y) : m_x(x), m_y(y) {}

/// @brief Create a copy of a 2D vector
/// @param other 2D vector
Vector2D::Vector2D(const Vector2D &other) : m_x(other.m_x), m_y(other.m_y) {}

/// @brief Get the i-th component of the vector
/// @param i 0-based index
/// @returns The i-th component
real_t Vector2D::operator[](const size_t i) const {
  switch (i) {
  case 0:
    return m_x;
  case 1:
    return m_y;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @brief Get the i-th component of the vector
/// @param i 0-based index
/// @returns The i-th component
real_t &Vector2D::operator[](const size_t i) {
  switch (i) {
  case 0:
    return m_x;
  case 1:
    return m_y;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @brief Get the i-th component of the vector
/// @param i 1-based index (may be positive or negative)
/// @returns The i-th component
real_t Vector2D::operator()(const index_t i) const {
  switch (i) {
  case 1:
    return m_x;
  case 2:
    return m_y;
  case -1:
    return m_y;
  case -2:
    return m_x;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @brief Get the i-th component of the vector
/// @param i 1-based index (may be positive or negative)
/// @returns The i-th component
real_t &Vector2D::operator()(const index_t i) {
  switch (i) {
  case 1:
    return m_x;
  case 2:
    return m_y;
  case -1:
    return m_y;
  case -2:
    return m_x;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @returns A copy of the 2D vector
Vector2D Vector2D::operator+() const { return Vector2D(*this); }

/// @brief Negate a 2D vector
/// @returns A negated copy of the 2D vector
Vector2D Vector2D::operator-() const {
  return Vector2D(-this->m_x, -this->m_y);
}

/// @brief Adds another 2D vector element-wise to this 2D vector
/// @param rhs 2D vector
/// @returns Updated 2D vector
Vector2D &Vector2D::operator+=(const Vector2D &rhs) {
  this->m_x += rhs.m_x;
  this->m_y += rhs.m_y;
  return *this;
}

/// @brief Subtracts another 2D vector element-wise from this 2D vector
/// @param rhs 2D vector
/// @returns Updated 2D vector
Vector2D &Vector2D::operator-=(const Vector2D &rhs) {
  this->m_x -= rhs.m_x;
  this->m_y -= rhs.m_y;
  return *this;
}

/// @brief Multiplies another 2D vector element-wise with this 2D vector
/// @param rhs 2D vector
/// @returns Updated 2D vector
Vector2D &Vector2D::operator*=(const Vector2D &rhs) {
  this->m_x *= rhs.m_x;
  this->m_y *= rhs.m_y;
  return *this;
}

/// @brief Divides this 2D vector by another 2D vector element-wise
/// @param rhs 2D vector
/// @returns Updated 2D vector
Vector2D &Vector2D::operator/=(const Vector2D &rhs) {
  if (isFuzzyEqual(rhs.m_x, 0.0) || isFuzzyEqual(rhs.m_y, 0.0)) {
    throw DivByZeroError();
  }
  this->m_x /= rhs.m_x;
  this->m_y /= rhs.m_y;
  return *this;
}

/// @brief Adds a scalar value element-wise to this 2D vector
/// @param rhs Scalar value
/// @returns Updated 2D vector
Vector2D &Vector2D::operator+=(real_t rhs) {
  this->m_x += rhs;
  this->m_y += rhs;
  return *this;
}

/// @brief Subtracts a scalar value element-wise from this 2D vector
/// @param rhs Scalar value
/// @returns Updated 2D vector
Vector2D &Vector2D::operator-=(real_t rhs) {
  this->m_x -= rhs;
  this->m_y -= rhs;
  return *this;
}

/// @brief Multiplies this 2D vector element-wise with a scalar value
/// @param rhs Scalar value
/// @returns Updated 2D vector
Vector2D &Vector2D::operator*=(real_t rhs) {
  this->m_x *= rhs;
  this->m_y *= rhs;
  return *this;
}

/// @brief Divides this vector element-wise by a scalar value
/// @param rhs Scalar value
/// @returns Updated 2D vector
Vector2D &Vector2D::operator/=(real_t rhs) {
  if (isFuzzyEqual(rhs, 0.0)) {
    throw DivByZeroError();
  }
  this->m_x /= rhs;
  this->m_y /= rhs;
  return *this;
}

/// @brief Get the x component of a 2D vector
/// @param vec 2D vector
/// @returns The x component
real_t x(const Vector2D &vec) { return vec.m_x; }

/// @brief Get the y component of a 2D vector
/// @param vec 2D vector
/// @returns The y component
real_t y(const Vector2D &vec) { return vec.m_y; }

/// @brief Computes the norm of a 2D vector
/// @param vec 2D Vector
/// @returns The vector's norm
real_t norm(const Vector2D &vec) { return std::sqrt(norm2(vec)); }

/// @brief Computes the squared norm of a 2D vector
/// @param vec 2D Vector
/// @returns The vector's squared norm
real_t norm2(const Vector2D &vec) { return sq(vec.m_x) + sq(vec.m_y); }

/// @brief Computes the absolute value of a 2D vector (element-wise)
/// @param vec 2D Vector
/// @returns The absolute value of a 2D vector (element-wise)
Vector2D abs(const Vector2D &vec) {
  return Vector2D(std::abs(vec.m_x), std::abs(vec.m_y));
}

/// @brief Computes the element-wise exponential of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise exponential
Vector2D exp(const Vector2D &vec) {
  return Vector2D(std::exp(vec.m_x), std::exp(vec.m_y));
}

/// @brief Computes the element-wise logarithm of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise logarithm
Vector2D log(const Vector2D &vec) {
  return Vector2D(std::log(vec.m_x), std::log(vec.m_y));
}

/// @brief Computes the element-wise sine of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise sine
Vector2D sin(const Vector2D &vec) {
  return Vector2D(std::sin(vec.m_x), std::sin(vec.m_y));
}

/// @brief Computes the element-wise cosine of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise cosine
Vector2D cos(const Vector2D &vec) {
  return Vector2D(std::cos(vec.m_x), std::cos(vec.m_y));
}

/// @brief Computes the element-wise tangent of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise tangent
Vector2D tan(const Vector2D &vec) {
  return Vector2D(std::tan(vec.m_x), std::tan(vec.m_y));
}

/// @brief Computes the element-wise arc-sine of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise arc-sine
Vector2D asin(const Vector2D &vec) {
  return Vector2D(std::asin(vec.m_x), std::asin(vec.m_y));
}

/// @brief Computes the element-wise arc-cosine of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise arc-cosine
Vector2D acos(const Vector2D &vec) {
  return Vector2D(std::acos(vec.m_x), std::acos(vec.m_y));
}

/// @brief Computes the element-wise arc-tangent of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise arc-tangent
Vector2D atan(const Vector2D &vec) {
  return Vector2D(std::atan(vec.m_x), std::atan(vec.m_y));
}

/// @brief Computes the element-wise power of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise power
Vector2D pow(const Vector2D &vec, real_t exponent) {
  return Vector2D(std::pow(vec.m_x, exponent), std::pow(vec.m_y, exponent));
}

/// @brief Computes the element-wise square root of a 2D vector
/// @param vec 2D vector
/// @returns The element-wise square root
Vector2D sqrt(const Vector2D &vec) {
  if (vec.m_x < 0.0 || vec.m_y < 0.0) {
    throw MathLibError("Square root of negative numbers not supported");
  }
  return Vector2D(std::sqrt(vec.m_x), std::sqrt(vec.m_y));
}

/// @brief Adds two 2D vectors element-wise
/// @param lhs 2D vector
/// @param rhs 2D vector
/// @returns Element-wise sum of the two vectors
Vector2D operator+(Vector2D lhs, const Vector2D &rhs) {
  lhs += rhs;
  return lhs;
}

/// @brief Adds a scalar scalar value element-wise to a 2D vector
/// @param lhs 2D vector
/// @param rhs Scalar value
/// @returns Element-wise sum of the 2D vector and sclar value
Vector2D operator+(Vector2D lhs, real_t rhs) {
  lhs += rhs;
  return lhs;
}

/// @brief Adds a scalar scalar value element-wise to a 2D vector
/// @param lhs Scalar value
/// @param rhs 2D vector
/// @returns Element-wise sum of the 2D vector and sclar value
Vector2D operator+(real_t lhs, Vector2D rhs) {
  rhs += lhs;
  return rhs;
}

/// @brief Subtracts two 2D vectors element-wise
/// @param lhs 2D vector
/// @param rhs 2D vector
/// @returns Element-wise subtraction of the two vectors
Vector2D operator-(Vector2D lhs, const Vector2D &rhs) {
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a scalar scalar value element-wise to a 2D vector
/// @param lhs 2D vector
/// @param rhs Scalar value
/// @returns Element-wise subtraction of the 2D vector and sclar value
Vector2D operator-(Vector2D lhs, real_t rhs) {
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a 2D vector from a scalar scalar value element-wise
/// @param lhs Scalar value
/// @param rhs 2D vector
/// @returns Element-wise subtraction of the sclar value and 2D vector
Vector2D operator-(real_t lhs, const Vector2D &rhs) {
  Vector2D tmp = -rhs;
  tmp += lhs;
  return tmp;
}

/// @brief Multiplies two 2D vectors element-wise
/// @param lhs 2D vector
/// @param rhs 2D vector
/// @returns Element-wise multiplication of the two vectors
real_t operator*(Vector2D lhs, const Vector2D &rhs) {
  return lhs.x() * rhs.x() + lhs.y() * rhs.y();
}

/// @brief Multiplies a scalar scalar value element-wise to a 2D vector
/// @param lhs 2D vector
/// @param rhs Scalar value
/// @returns Element-wise multiplication of the 2D vector and sclar value
Vector2D operator*(Vector2D lhs, real_t rhs) {
  lhs *= rhs;
  return lhs;
}

/// @brief Multiplies a scalar scalar value element-wise to a 2D vector
/// @param lhs Scalar value
/// @param rhs 2D vector
/// @returns Element-wise multiplication of the 2D vector and sclar value
Vector2D operator*(real_t lhs, Vector2D rhs) {
  rhs *= lhs;
  return rhs;
}

/// @brief Divides a 2D vector element-wise by a scalar value
/// @param lhs 2D vector
/// @param rhs Scalar value
/// @returns Element-wise division of the 2D vector by the scalar value
Vector2D operator/(Vector2D lhs, real_t rhs) {
  if (isFuzzyEqual(rhs, 0.0)) {
    throw DivByZeroError();
  }
  lhs /= rhs;
  return lhs;
}

/// @brief Divides a scalar value element-wise by a 2D vector
/// @param lhs Scalar value
/// @param rhs 2D vector
/// @returns Element-wise division of the scalar value by the 2D vector
Vector2D operator/(real_t lhs, const Vector2D &rhs) {
  if (isFuzzyEqual(rhs.x(), 0.0) || isFuzzyEqual(rhs.y(), 0.0)) {
    throw DivByZeroError();
  }
  return Vector2D(lhs / rhs.x(), lhs / rhs.y());
}

/// @brief Compute the dot product of two 2D vectors
/// @param vec1 2D vector
/// @param vec2 2D vector
/// @return Dot product of the two vectors
real_t dot(const Vector2D &vec1, const Vector2D &vec2) { return vec1 * vec2; }

/// @brief Normalize a vector
/// @param vec 2D vector
/// @return A normalized copy of the vector
Vector2D normalize(const Vector2D &vec) { return vec / norm(vec); }

/// @brief Compute the angle between a vector and the x axis
/// @param vec 2D vector
/// @return Angle [rad]
real_t angle(const Vector2D &vec) { return std::atan2(vec.y(), vec.x()); }

/// @brief Compute the angle between two vectors
/// @param vec1 2D vector
/// @param vec2 2D vector
/// @return Angle [rad]
real_t angle(const Vector2D &vec1, const Vector2D &vec2) {
  return angle(vec2) - angle(vec1);
}

#ifdef __PRINT__

/// @brief Converts a 2D vector to a string
/// @param vec 2D Vector
/// @returns The complex number as a string
std::string to_string(const Vector2D &vec, const real_t epsilon) {
  std::string str = "[";
  str += std::to_string(vec.x());
  str += ", ";
  str += std::to_string(vec.y());
  str += "]";
  return str;
}

/// @brief Serializes a 2D vector
std::ostream &operator<<(std::ostream &os, const Vector2D &vec) {
  os << to_string(vec);
  return os;
}

/// @brief Prints a 2D vector in the default output stream
void print(const Vector2D &vec) { std::cout << vec << '\n'; }

#endif

} // namespace mathlib

#ifdef __PRINT__
namespace std {

/// @brief Overloads the conversion from a 2D vector to a string
std::string to_string(const mathlib::Vector2D &vec) {
  return mathlib::to_string(vec);
}
#endif

} // namespace std