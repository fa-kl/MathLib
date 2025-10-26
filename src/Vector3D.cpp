/*****************************************************************************************
 * @file: Vector3D.hpp
 *
 * @brief: A class for 3D vectors.
 *
 * @details: This file implements a class Vector3D for handling 3D vectors.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "Vector3D.hpp"
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

/// @brief Create a 3D vector of all zeros
Vector3D::Vector3D() : m_x(0.0), m_y(0.0), m_z(0.0) {}

/// @brief Create a 3D vector
/// @param x The vector's x component
/// @param y The vector's y component
/// @param z The vector's z component
Vector3D::Vector3D(const real_t x, const real_t y, const real_t z)
    : m_x(x), m_y(y), m_z(z) {}

/// @brief Create a copy of a 3D vector
/// @param other 3D vector
Vector3D::Vector3D(const Vector3D &other)
    : m_x(other.m_x), m_y(other.m_y), m_z(other.m_z) {}

/// @brief Get the i-th component of the vector
/// @param i 0-based index
/// @returns The i-th component
real_t Vector3D::operator[](const size_t i) const {
  switch (i) {
  case 0:
    return m_x;
  case 1:
    return m_y;
  case 2:
    return m_z;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @brief Get the i-th component of the vector
/// @param i 0-based index
/// @returns The i-th component
real_t &Vector3D::operator[](const size_t i) {
  switch (i) {
  case 0:
    return m_x;
  case 1:
    return m_y;
  case 2:
    return m_z;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @brief Get the i-th component of the vector
/// @param i 1-based index (may be positive or negative)
/// @returns The i-th component
real_t Vector3D::operator()(const index_t i) const {
  switch (i) {
  case 1:
    return m_x;
  case 2:
    return m_y;
  case 3:
    return m_z;
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
real_t &Vector3D::operator()(const index_t i) {
  switch (i) {
  case 1:
    return m_x;
  case 2:
    return m_y;
  case 3:
    return m_z;
  case -1:
    return m_y;
  case -2:
    return m_x;
  default:
    throw IndexOutOfRangeError();
  }
  throw MathLibError("Well that's unexpected.");
}

/// @returns A copy of the 3D vector
Vector3D Vector3D::operator+() const { return Vector3D(*this); }

/// @brief Negate a 3D vector
/// @returns A negated copy of the 3D vector
Vector3D Vector3D::operator-() const {
  return Vector3D(-this->m_x, -this->m_y, -this->m_z);
}

/// @brief Adds another 3D vector element-wise to this 3D vector
/// @param rhs 3D vector
/// @returns Updated 3D vector
Vector3D &Vector3D::operator+=(const Vector3D &rhs) {
  this->m_x += rhs.m_x;
  this->m_y += rhs.m_y;
  this->m_z += rhs.m_z;
  return *this;
}

/// @brief Subtracts another 3D vector element-wise from this 3D vector
/// @param rhs 3D vector
/// @returns Updated 3D vector
Vector3D &Vector3D::operator-=(const Vector3D &rhs) {
  this->m_x -= rhs.m_x;
  this->m_y -= rhs.m_y;
  this->m_z -= rhs.m_z;
  return *this;
}

/// @brief Multiplies another 3D vector element-wise with this 3D vector
/// @param rhs 3D vector
/// @returns Updated 3D vector
Vector3D &Vector3D::operator*=(const Vector3D &rhs) {
  this->m_x *= rhs.m_x;
  this->m_y *= rhs.m_y;
  this->m_z *= rhs.m_z;
  return *this;
}

/// @brief Divides this 3D vector by another 3D vector element-wise
/// @param rhs 3D vector
/// @returns Updated 3D vector
Vector3D &Vector3D::operator/=(const Vector3D &rhs) {
  if (isFuzzyEqual(rhs.m_x, 0.0) || isFuzzyEqual(rhs.m_y, 0.0) ||
      isFuzzyEqual(rhs.m_z, 0.0)) {
    throw DivByZeroError();
  }
  this->m_x /= rhs.m_x;
  this->m_y /= rhs.m_y;
  this->m_z /= rhs.m_z;
  return *this;
}

/// @brief Adds a scalar value element-wise to this 3D vector
/// @param rhs Scalar value
/// @returns Updated 3D vector
Vector3D &Vector3D::operator+=(real_t rhs) {
  this->m_x += rhs;
  this->m_y += rhs;
  this->m_z += rhs;
  return *this;
}

/// @brief Subtracts a scalar value element-wise from this 3D vector
/// @param rhs Scalar value
/// @returns Updated 3D vector
Vector3D &Vector3D::operator-=(real_t rhs) {
  this->m_x -= rhs;
  this->m_y -= rhs;
  this->m_z -= rhs;
  return *this;
}

/// @brief Multiplies this 3D vector element-wise with a scalar value
/// @param rhs Scalar value
/// @returns Updated 3D vector
Vector3D &Vector3D::operator*=(real_t rhs) {
  this->m_x *= rhs;
  this->m_y *= rhs;
  this->m_z *= rhs;
  return *this;
}

/// @brief Divides this vector element-wise by a scalar value
/// @param rhs Scalar value
/// @returns Updated 3D vector
Vector3D &Vector3D::operator/=(real_t rhs) {
  if (isFuzzyEqual(rhs, 0.0)) {
    throw DivByZeroError();
  }
  this->m_x /= rhs;
  this->m_y /= rhs;
  this->m_z /= rhs;
  return *this;
}

/// @brief Get the x component of a 3D vector
/// @param vec 3D vector
/// @returns The x component
real_t x(const Vector3D &vec) { return vec.m_x; }

/// @brief Get the y component of a 3D vector
/// @param vec 3D vector
/// @returns The y component
real_t y(const Vector3D &vec) { return vec.m_y; }

/// @brief Get the z component of a 3D vector
/// @param vec 3D vector
/// @returns The z component
real_t z(const Vector3D &vec) { return vec.m_z; }

/// @brief Computes the norm of a 3D vector
/// @param vec 3D Vector
/// @returns The vector's norm
real_t norm(const Vector3D &vec) { return std::sqrt(norm2(vec)); }

/// @brief Computes the squared norm of a 3D vector
/// @param vec 3D Vector
/// @returns The vector's squared norm
real_t norm2(const Vector3D &vec) {
  return sq(vec.m_x) + sq(vec.m_y) + sq(vec.m_z);
}

/// @brief Computes the absolute value of a 3D vector (element-wise)
/// @param vec 3D Vector
/// @returns The absolute value of a 3D vector (element-wise)
Vector3D abs(const Vector3D &vec) {
  return Vector3D(std::abs(vec.m_x), std::abs(vec.m_y), std::abs(vec.m_z));
}

/// @brief Computes the element-wise exponential of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise exponential
Vector3D exp(const Vector3D &vec) {
  return Vector3D(std::exp(vec.m_x), std::exp(vec.m_y), std::exp(vec.m_z));
}

/// @brief Computes the element-wise logarithm of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise logarithm
Vector3D log(const Vector3D &vec) {
  return Vector3D(std::log(vec.m_x), std::log(vec.m_y), std::log(vec.m_z));
}

/// @brief Computes the element-wise sine of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise sine
Vector3D sin(const Vector3D &vec) {
  return Vector3D(std::sin(vec.m_x), std::sin(vec.m_y), std::sin(vec.m_z));
}

/// @brief Computes the element-wise cosine of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise cosine
Vector3D cos(const Vector3D &vec) {
  return Vector3D(std::cos(vec.m_x), std::cos(vec.m_y), std::cos(vec.m_z));
}

/// @brief Computes the element-wise tangent of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise tangent
Vector3D tan(const Vector3D &vec) {
  return Vector3D(std::tan(vec.m_x), std::tan(vec.m_y), std::tan(vec.m_z));
}

/// @brief Computes the element-wise arc-sine of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise arc-sine
Vector3D asin(const Vector3D &vec) {
  return Vector3D(std::asin(vec.m_x), std::asin(vec.m_y), std::asin(vec.m_z));
}

/// @brief Computes the element-wise arc-cosine of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise arc-cosine
Vector3D acos(const Vector3D &vec) {
  return Vector3D(std::acos(vec.m_x), std::acos(vec.m_y), std::acos(vec.m_z));
}

/// @brief Computes the element-wise arc-tangent of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise arc-tangent
Vector3D atan(const Vector3D &vec) {
  return Vector3D(std::atan(vec.m_x), std::atan(vec.m_y), std::atan(vec.m_z));
}

/// @brief Computes the element-wise power of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise power
Vector3D pow(const Vector3D &vec, real_t exponent) {
  return Vector3D(std::pow(vec.m_x, exponent), std::pow(vec.m_y, exponent),
                  std::pow(vec.m_z, exponent));
}

/// @brief Computes the element-wise square root of a 3D vector
/// @param vec 3D vector
/// @returns The element-wise square root
Vector3D sqrt(const Vector3D &vec) {
  if (vec.m_x < 0.0 || vec.m_y < 0.0 || vec.m_z < 0.0) {
    throw MathLibError("Square root of negative numbers not supported");
  }
  return Vector3D(std::sqrt(vec.m_x), std::sqrt(vec.m_y), std::sqrt(vec.m_z));
}

/// @brief Adds two 3D vectors element-wise
/// @param lhs 3D vector
/// @param rhs 3D vector
/// @returns Element-wise sum of the two vectors
Vector3D operator+(Vector3D lhs, const Vector3D &rhs) {
  lhs += rhs;
  return lhs;
}

/// @brief Adds a scalar scalar value element-wise to a 3D vector
/// @param lhs 3D vector
/// @param rhs Scalar value
/// @returns Element-wise sum of the 3D vector and sclar value
Vector3D operator+(Vector3D lhs, real_t rhs) {
  lhs += rhs;
  return lhs;
}

/// @brief Adds a scalar scalar value element-wise to a 3D vector
/// @param lhs Scalar value
/// @param rhs 3D vector
/// @returns Element-wise sum of the 3D vector and sclar value
Vector3D operator+(real_t lhs, Vector3D rhs) {
  rhs += lhs;
  return rhs;
}

/// @brief Subtracts two 3D vectors element-wise
/// @param lhs 3D vector
/// @param rhs 3D vector
/// @returns Element-wise subtraction of the two vectors
Vector3D operator-(Vector3D lhs, const Vector3D &rhs) {
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a scalar scalar value element-wise to a 3D vector
/// @param lhs 3D vector
/// @param rhs Scalar value
/// @returns Element-wise subtraction of the 3D vector and sclar value
Vector3D operator-(Vector3D lhs, real_t rhs) {
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a 3D vector from a scalar scalar value element-wise
/// @param lhs Scalar value
/// @param rhs 3D vector
/// @returns Element-wise subtraction of the sclar value and 3D vector
Vector3D operator-(real_t lhs, const Vector3D &rhs) {
  Vector3D tmp = -rhs;
  tmp += lhs;
  return tmp;
}

/// @brief Compute the dot product of two 3D vectors
/// @param vec1 3D vector
/// @param vec2 3D vector
/// @return Dot product of the two vectors
real_t operator*(Vector3D lhs, const Vector3D &rhs) {
  return lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z();
}

/// @brief Multiplies a scalar scalar value element-wise to a 3D vector
/// @param lhs 3D vector
/// @param rhs Scalar value
/// @returns Element-wise multiplication of the 3D vector and sclar value
Vector3D operator*(Vector3D lhs, real_t rhs) {
  lhs *= rhs;
  return lhs;
}

/// @brief Multiplies a scalar scalar value element-wise to a 3D vector
/// @param lhs Scalar value
/// @param rhs 3D vector
/// @returns Element-wise multiplication of the 3D vector and sclar value
Vector3D operator*(real_t lhs, Vector3D rhs) {
  rhs *= lhs;
  return rhs;
}

/// @brief Divides a 3D vector element-wise by a scalar value
/// @param lhs 3D vector
/// @param rhs Scalar value
/// @returns Element-wise division of the 3D vector by the scalar value
Vector3D operator/(Vector3D lhs, real_t rhs) {
  if (isFuzzyEqual(rhs, 0.0)) {
    throw DivByZeroError();
  }
  lhs /= rhs;
  return lhs;
}

/// @brief Divides a scalar value element-wise by a 3D vector
/// @param lhs Scalar value
/// @param rhs 3D vector
/// @returns Element-wise division of the scalar value by the 3D vector
Vector3D operator/(real_t lhs, const Vector3D &rhs) {
  if (isFuzzyEqual(rhs.x(), 0.0) || isFuzzyEqual(rhs.y(), 0.0) ||
      isFuzzyEqual(rhs.z(), 0.0)) {
    throw DivByZeroError();
  }
  return Vector3D(lhs / rhs.x(), lhs / rhs.y(), lhs / rhs.z());
}

/// @brief Compute the dot product of two 3D vectors
/// @param vec1 3D vector
/// @param vec2 3D vector
/// @return Dot product of the two vectors
real_t dot(const Vector3D &vec1, const Vector3D &vec2) { return vec1 * vec2; }

/// @brief Compute the cross product of two 3D vectors
/// @param vec1 3D vector
/// @param vec2 3D vector
/// @return Cross product of the two vectors
Vector3D cross(const Vector3D &vec1, const Vector3D &vec2) {
  return Vector3D(vec1.y() * vec2.z() - vec1.z() * vec2.y(),
                  vec1.z() * vec2.x() - vec1.x() * vec2.z(),
                  vec1.x() * vec2.y() - vec1.y() * vec2.x());
}

/// @brief Normalize a vector
/// @param vec 3D vector
/// @return A normalized copy of the vector
Vector3D normalize(const Vector3D &vec) { return vec / norm(vec); }

#ifdef __PRINT__

/// @brief Converts a 3D vector to a string
/// @param vec 3D Vector
/// @returns The complex number as a string
std::string to_string(const Vector3D &vec, const real_t epsilon) {
  std::string str = "[";
  str += std::to_string(vec.x());
  str += ", ";
  str += std::to_string(vec.y());
  str += ", ";
  str += std::to_string(vec.z());
  str += "]";
  return str;
}

/// @brief Serializes a 3D vector
std::ostream &operator<<(std::ostream &os, const Vector3D &vec) {
  os << to_string(vec);
  return os;
}

/// @brief Prints a 3D vector in the default output stream
void print(const Vector3D &vec) { std::cout << vec << '\n'; }

#endif

} // namespace mathlib

#ifdef __PRINT__
namespace std {

/// @brief Overloads the conversion from a 3D vector to a string
std::string to_string(const mathlib::Vector3D &vec) {
  return mathlib::to_string(vec);
}
#endif

} // namespace std