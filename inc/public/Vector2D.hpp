/*****************************************************************************************
 * @file: Vector2D.hpp
 *
 * @brief: A class for 2D vectors.
 *
 * @details: This file declares a class Vector2D for handling 2D vectors.
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

/// @brief A class for 2D vectors
class Vector2D {
protected:
  /// @brief The x component of the 2D vector
  real_t m_x;
  /// @brief The y component of the 2D vector
  real_t m_y;

public:
  Vector2D();
  Vector2D(real_t x, real_t y);
  Vector2D(const Vector2D &other);

  Vector2D &operator=(const Vector2D &other) = default;

  /// Get the x component of the vector
  /// @returns The x component
  real_t x() const { return m_x; }

  /// Get the y component of the vector
  /// @returns The y component
  real_t y() const { return m_y; }

  /// Get the x component of the vector
  /// @returns The x component
  real_t &x() { return m_x; }

  /// Get the y component of the vector
  /// @returns The y component
  real_t &y() { return m_y; }

  real_t operator[](const size_t i) const;
  real_t &operator[](const size_t i);
  real_t operator()(const index_t i) const;
  real_t &operator()(const index_t i);

  Vector2D operator+() const;
  Vector2D operator-() const;
  Vector2D &operator+=(const Vector2D &rhs);
  Vector2D &operator-=(const Vector2D &rhs);
  Vector2D &operator*=(const Vector2D &rhs);
  Vector2D &operator/=(const Vector2D &rhs);
  Vector2D &operator+=(real_t rhs);
  Vector2D &operator-=(real_t rhs);
  Vector2D &operator*=(real_t rhs);
  Vector2D &operator/=(real_t rhs);

  friend real_t x(const Vector2D &vec);
  friend real_t y(const Vector2D &vec);
  friend real_t norm(const Vector2D &vec);
  friend real_t norm2(const Vector2D &vec);
  friend Vector2D abs(const Vector2D &vec);
  friend Vector2D exp(const Vector2D &vec);
  friend Vector2D log(const Vector2D &vec);
  friend Vector2D sin(const Vector2D &vec);
  friend Vector2D cos(const Vector2D &vec);
  friend Vector2D tan(const Vector2D &vec);
  friend Vector2D asin(const Vector2D &vec);
  friend Vector2D acos(const Vector2D &vec);
  friend Vector2D atan(const Vector2D &vec);
  friend Vector2D pow(const Vector2D &vec, real_t exponent);
  friend Vector2D sqrt(const Vector2D &vec);
  friend real_t dot(const Vector2D &vec1, const Vector2D &vec2);
  friend Vector2D normalize(const Vector2D &vec);
  friend real_t angle(const Vector2D &vec);
  friend real_t angle(const Vector2D &vec1, const Vector2D &vec2);
};

Vector2D operator+(Vector2D lhs, const Vector2D &rhs);
Vector2D operator+(Vector2D lhs, real_t rhs);
Vector2D operator+(real_t lhs, Vector2D rhs);
Vector2D operator-(Vector2D lhs, const Vector2D &rhs);
Vector2D operator-(Vector2D lhs, real_t rhs);
Vector2D operator-(real_t lhs, const Vector2D &rhs);
real_t operator*(Vector2D lhs, const Vector2D &rhs);
Vector2D operator*(Vector2D lhs, real_t rhs);
Vector2D operator*(real_t lhs, Vector2D rhs);
Vector2D operator/(Vector2D lhs, real_t rhs);
Vector2D operator/(real_t lhs, const Vector2D &rhs);

#ifdef __PRINT__
std::string to_string(const Vector2D &vec, const real_t epsilon = 1e-6);
std::ostream &operator<<(std::ostream &os, const Vector2D &vec);
void print(const Vector2D &vec);
#endif

} // namespace mathlib

namespace std {

std::string to_string(const mathlib::Vector2D &vec);

} // namespace std
