/*****************************************************************************************
 * @file: Complex.hpp
 *
 * @brief: A class for complex numbers.
 *
 * @details: This file declares a class Complex for handling complex numbers.
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

/// @brief A class for complex numbers
class Complex {
protected:
  /// @brief The real component of the complex number
  real_t m_real;
  /// @brief The imaginary component of the complex number
  real_t m_imag;

public:
  Complex();
  Complex(real_t real, real_t imag);
  Complex(const Complex &other);

  Complex &operator=(const Complex &other);

  /// Get the real part of the complex number
  /// @returns Real part
  real_t real() const { return m_real; }

  /// Get the imaginary part of the complex number
  /// @returns Imaginary part
  real_t imag() const { return m_imag; }

  Complex operator+() const;
  Complex operator-() const;
  Complex &operator+=(const Complex &rhs);
  Complex &operator-=(const Complex &rhs);
  Complex &operator*=(const Complex &rhs);
  Complex &operator/=(const Complex &rhs);
  Complex &operator+=(real_t rhs);
  Complex &operator-=(real_t rhs);
  Complex &operator*=(real_t rhs);
  Complex &operator/=(real_t rhs);

  friend real_t real(const Complex &z);
  friend real_t imag(const Complex &z);
  friend real_t abs(const Complex &z);
  friend real_t abs2(const Complex &z);
  friend real_t arg(const Complex &z);
  friend Complex conj(const Complex &z);
  friend Complex exp(const Complex &z);
  friend Complex log(const Complex &z);
  friend Complex sin(const Complex &z);
  friend Complex cos(const Complex &z);
  friend Complex tan(const Complex &z);
  friend Complex asin(const Complex &z);
  friend Complex acos(const Complex &z);
  friend Complex atan(const Complex &z);
  friend Complex pow(const Complex &z, real_t w);
  friend Complex pow(const Complex &z, const Complex &w);
  friend Complex sqrt(const Complex &z);
};

Complex operator+(Complex lhs, const Complex &rhs);
Complex operator+(Complex lhs, real_t rhs);
Complex operator+(real_t lhs, Complex rhs);
Complex operator-(Complex lhs, const Complex &rhs);
Complex operator-(Complex lhs, real_t rhs);
Complex operator-(real_t lhs, const Complex &rhs);
Complex operator*(Complex lhs, const Complex &rhs);
Complex operator*(Complex lhs, real_t rhs);
Complex operator*(real_t lhs, Complex rhs);
Complex operator/(Complex lhs, const Complex &rhs);
Complex operator/(Complex lhs, real_t rhs);
Complex operator/(real_t lhs, const Complex &rhs);

#ifdef __PRINT__
std::string to_string(const Complex &z, const real_t epsilon = 1e-6);
std::ostream &operator<<(std::ostream &os, const Complex &z);
void print(const Complex &z);
#endif

/// @brief The imaginary unit
extern const Complex j;

} // namespace mathlib

namespace std {

std::string to_string(const mathlib::Complex &z);

} // namespace std
