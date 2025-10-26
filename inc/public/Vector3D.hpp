/*****************************************************************************************
 * @file: Vector3D.hpp
 *
 * @brief: A class for 3D vectors.
 *
 * @details: This file declares a class Vector3D for handling 3D vectors.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include "config.hpp"

#ifdef __PRINT__
#include <iosfwd>
#include <string>
#endif

namespace mathlib {

/// @brief A class for 3D vectors
class Vector3D {
protected:
  /// @brief The x component of the 3D vector
  real_t m_x;
  /// @brief The y component of the 3D vector
  real_t m_y;
  /// @brief The z component of the 3D vector
  real_t m_z;

public:
  Vector3D();
  Vector3D(const real_t x, const real_t y, const real_t z);
  Vector3D(const Vector3D &other);

  Vector3D &operator=(const Vector3D &other) = default;

  /// Get the x component of the vector
  /// @returns The x component
  real_t x() const { return m_x; }

  /// Get the y component of the vector
  /// @returns The y component
  real_t y() const { return m_y; }

  /// Get the z component of the vector
  /// @returns The z component
  real_t z() const { return m_z; }

  /// Get the x component of the vector
  /// @returns The x component
  real_t &x() { return m_x; }

  /// Get the y component of the vector
  /// @returns The y component
  real_t &y() { return m_y; }

  /// Get the z component of the vector
  /// @returns The z component
  real_t &z() { return m_z; }

  real_t operator[](const size_t i) const;
  real_t &operator[](const size_t i);
  real_t operator()(const index_t i) const;
  real_t &operator()(const index_t i);

  Vector3D operator+() const;
  Vector3D operator-() const;
  Vector3D &operator+=(const Vector3D &rhs);
  Vector3D &operator-=(const Vector3D &rhs);
  Vector3D &operator*=(const Vector3D &rhs);
  Vector3D &operator/=(const Vector3D &rhs);
  Vector3D &operator+=(real_t rhs);
  Vector3D &operator-=(real_t rhs);
  Vector3D &operator*=(real_t rhs);
  Vector3D &operator/=(real_t rhs);

  friend real_t x(const Vector3D &vec);
  friend real_t y(const Vector3D &vec);
  friend real_t z(const Vector3D &vec);
  friend real_t norm(const Vector3D &vec);
  friend real_t norm2(const Vector3D &vec);
  friend Vector3D abs(const Vector3D &vec);
  friend Vector3D exp(const Vector3D &vec);
  friend Vector3D log(const Vector3D &vec);
  friend Vector3D sin(const Vector3D &vec);
  friend Vector3D cos(const Vector3D &vec);
  friend Vector3D tan(const Vector3D &vec);
  friend Vector3D asin(const Vector3D &vec);
  friend Vector3D acos(const Vector3D &vec);
  friend Vector3D atan(const Vector3D &vec);
  friend Vector3D pow(const Vector3D &vec, real_t exponent);
  friend Vector3D sqrt(const Vector3D &vec);
  friend real_t dot(const Vector3D &vec1, const Vector3D &vec2);
  friend Vector3D cross(const Vector3D &vec1, const Vector3D &vec2);
  friend Vector3D normalize(const Vector3D &vec);
};

Vector3D operator+(Vector3D lhs, const Vector3D &rhs);
Vector3D operator+(Vector3D lhs, real_t rhs);
Vector3D operator+(real_t lhs, Vector3D rhs);
Vector3D operator-(Vector3D lhs, const Vector3D &rhs);
Vector3D operator-(Vector3D lhs, real_t rhs);
Vector3D operator-(real_t lhs, const Vector3D &rhs);
real_t operator*(Vector3D lhs, const Vector3D &rhs);
Vector3D operator*(Vector3D lhs, real_t rhs);
Vector3D operator*(real_t lhs, Vector3D rhs);
Vector3D operator/(Vector3D lhs, real_t rhs);
Vector3D operator/(real_t lhs, const Vector3D &rhs);

#ifdef __PRINT__
std::string to_string(const Vector3D &vec, const real_t epsilon = 1e-6);
std::ostream &operator<<(std::ostream &os, const Vector3D &vec);
void print(const Vector3D &vec);
#endif

} // namespace mathlib

namespace std {

std::string to_string(const mathlib::Vector3D &vec);

} // namespace std
