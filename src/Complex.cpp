/*****************************************************************************************
 * @file: Complex.cpp
 *
 * @brief: A class for complex numbers.
 *
 * @details: This file implements a class Complex for handling complex numbers.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "Complex.hpp"

#include <cmath>
#ifdef __PRINT__
#include <iostream>
#endif

namespace mathlib
{

const Complex j(0.0, 1.0);

/// @brief Creates a complex number with 0 real and 0 imaginary parts
Complex::Complex() : m_real(0.0), m_imag(0.0) {}

/// @brief Create a new complex number
/// @param real Real part
/// @param imag Imaginary part
Complex::Complex(const real_t real, const real_t imag) : m_real(real), m_imag(imag) {}

/// @brief Create a copy of a complex number
/// @param other Complex number
Complex::Complex(const Complex& other) : m_real(other.m_real), m_imag(other.m_imag) {}

/// @brief Copy operator
/// @param other Complex number
/// @returns Copied complex number
Complex& Complex::operator=(const Complex& other)
{
  if (this != &other) {
    m_real = other.m_real;
    m_imag = other.m_imag;
  }
  return *this;
}

/// @brief Get the real part of the complex number z
/// @param z Complex number
/// @returns Real part of z
real_t real(const Complex& z)
{
  return z.m_real;
}

/// @brief Get the imaginary part of the complex number z
/// @param z Complex number
/// @returns Imaginary part of z
real_t imag(const Complex& z)
{
  return z.m_imag;
}

/// @brief Computes the magnitude of a complex number z
/// @param z Complex number
/// @returns Magnitude of z
real_t abs(const Complex& z)
{
  return std::sqrt(abs2(z));
}

/// @brief Computes the squared magnitude of a complex number z
/// @param z Complex number
/// @returns Squared magnitude of z
real_t abs2(const Complex& z)
{
  return z.m_real * z.m_real + z.m_imag * z.m_imag;
}

/// @brief Computes the argument (angle to real axis) of a complex number z
/// @param z Complex number
/// @returns Argument of z
real_t arg(const Complex& z)
{
  return std::atan2(z.m_imag, z.m_real);
}

/// @brief Computes the conjugate of a complex number z
/// @param z Complex number
/// @returns Conjugate complex number of z
Complex conj(const Complex& z)
{
  return Complex(z.m_real, -z.m_imag);
}

/// @brief Computes the exponential map of a complex number z
/// @param z Complex number
/// @returns Exponential of z
Complex exp(const Complex& z)
{
  return std::exp(z.m_real) * (std::cos(z.m_imag) + j * std::sin(z.m_imag));
}

/// @brief Computes the logarithmic map of a complex number z
/// @param z Complex number
/// @returns Logarithm of z
Complex log(const Complex& z)
{
  return Complex(std::log(abs(z)), arg(z));
}

/// @brief Computes the sine of a complex number z
/// @param z Complex number
/// @returns Sine of z
Complex sin(const Complex& z)
{
  return Complex(std::sin(z.m_real) * std::cosh(z.m_imag), std::cos(z.m_real) * std::sinh(z.m_imag));
}

/// @brief Computes the cosine of a complex number z
/// @param z Complex number
/// @returns Cosine of z
Complex cos(const Complex& z)
{
  return Complex(std::cos(z.m_real) * std::cosh(z.m_imag), -std::sin(z.m_real) * std::sinh(z.m_imag));
}

/// @brief Computes the tangent of a complex number z
/// @param z Complex number
/// @returns Tangent of z
Complex tan(const Complex& z)
{
  return sin(z) / cos(z);
}

/// @brief Computes the arc-sine of a complex number z
/// @param z Complex number
/// @returns Arc-sine of z
Complex asin(const Complex& z)
{
  return -j * log(j * z + sqrt(static_cast<real_t>(1.0) - z * z));
}

/// @brief Computes the arc-cosine of a complex number z
/// @param z Complex number
/// @returns Arc-cosine of z
Complex acos(const Complex& z)
{
  return -j * log(z + sqrt(z * z - static_cast<real_t>(1.0)));
}

/// @brief Computes the arc-tangent of a complex number z
/// @param z Complex number
/// @returns Arc-tangent of z
Complex atan(const Complex& z)
{
  const Complex one(1, 0);
  return (j / Complex(2, 0)) * (log(one - j * z) - log(one + j * z));
}

/// @brief Computes the power of a complex number z with a real exponent w
/// @param z Base, complex number
/// @param w Exponent, real number
/// @returns z to the power of w
Complex pow(const Complex& z, const real_t w)
{
  return exp(w * log(z));
}

/// @brief Computes the power of a complex number z with a complex exponent w
/// @param z Base, complex number
/// @param w Exponent, complex number
/// @returns z to the power of w
Complex pow(const Complex& z, const Complex& w)
{
  return exp(w * log(z));
}

/// @brief Computes the square root of a complex number z
/// @param z Complex number
/// @returns Square root of z
Complex sqrt(const Complex& z)
{
  const real_t r = abs(z);
  const real_t x = z.m_real;
  const real_t y = z.m_imag;
  const real_t u = std::sqrt((r + x) / 2);
  real_t v = std::sqrt((r - x) / 2);
  if (y < 0) {
    v = -v;
  }
  return Complex(u, v);
}

/// @returns A copy of the complex number
Complex Complex::operator+() const
{
  return Complex(*this);
}

/// @brief Negate a complex number
/// @returns A negated copy
Complex Complex::operator-() const
{
  return Complex(-this->m_real, -this->m_imag);
}

/// @brief Adds two complex numbers and assigns it to the original variable
/// @param rhs Complex number
/// @returns The updated complex number
Complex& Complex::operator+=(const Complex& rhs)
{
  this->m_real += rhs.m_real;
  this->m_imag += rhs.m_imag;
  return *this;
}

/// @brief Subtracts two complex numbers and assigns it to the original variable
/// @param rhs The complex number to subtract
/// @returns The updated complex number
Complex& Complex::operator-=(const Complex& rhs)
{
  this->m_real -= rhs.m_real;
  this->m_imag -= rhs.m_imag;
  return *this;
}

/// @brief Multiplies two complex numbers and assigns it to the original
/// variable
/// @param rhs Complex number
/// @returns The updated complex number
Complex& Complex::operator*=(const Complex& rhs)
{
  const real_t real = this->m_real * rhs.m_real - this->m_imag * rhs.m_imag;
  this->m_imag = this->m_imag * rhs.m_real + this->m_real * rhs.m_imag;
  this->m_real = real;
  return *this;
}

/// @brief Divides two complex numbers and assigns it to the original variable
/// @param rhs Divisor, complex number
/// @returns The updated complex number
Complex& Complex::operator/=(const Complex& rhs)
{
  const real_t div = abs2(rhs);
  const real_t real = (this->m_real * rhs.m_real + this->m_imag * rhs.m_imag) / div;
  this->m_imag = (this->m_imag * rhs.m_real - this->m_real * rhs.m_imag) / div;
  this->m_real = real;
  return *this;
}

/// @brief Adds a real number to a complex number and assigns it to the original
/// variable
/// @param rhs Real number
/// @returns The updated complex number
Complex& Complex::operator+=(const real_t rhs)
{
  this->m_real += rhs;
  return *this;
}

/// @brief Subtracts a real number from a complex number and assigns it to the
/// original variable
/// @param rhs Real number
/// @returns The updated complex number
Complex& Complex::operator-=(const real_t rhs)
{
  this->m_real -= rhs;
  return *this;
}

/// @brief Multiplies a complex number with a real number and assigns it to the
/// original variable
/// @param rhs Real number
/// @returns The updated complex number
Complex& Complex::operator*=(const real_t rhs)
{
  this->m_real *= rhs;
  this->m_imag *= rhs;
  return *this;
}

/// @brief Divides a complex number by a real number  and assigns it to the
/// original variable
/// @param rhs Real number
/// @returns The updated complex number
Complex& Complex::operator/=(const real_t rhs)
{
  this->m_real /= rhs;
  this->m_imag /= rhs;
  return *this;
}

/// @brief Adds two complex numbers
/// @param lhs Complex number
/// @param rhs Complex number
/// @returns The addition of the two complex numbers
Complex operator+(Complex lhs, const Complex& rhs)
{
  lhs += rhs;
  return lhs;
}

/// @brief Adds a real number to a complex number
/// @param lhs Complex number
/// @param rhs Real Number
/// @returns The addition of the complex and real numbers
Complex operator+(Complex lhs, const real_t rhs)
{
  lhs += rhs;
  return lhs;
}

/// @brief Adds a real number to a complex number
/// @param lhs Real Number
/// @param rhs Complex number
/// @returns The addition of the complex and real numbers
Complex operator+(const real_t lhs, Complex rhs)
{
  rhs += lhs;
  return rhs;
}

/// @brief Subtracts two complex numbers
/// @param lhs Complex number
/// @param rhs Complex number
/// @returns The difference of the two complex numbers
Complex operator-(Complex lhs, const Complex& rhs)
{
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a real number from a complex number
/// @param lhs Complex number
/// @param rhs Real number
/// @returns The difference of the two complex numbers
Complex operator-(Complex lhs, const real_t rhs)
{
  lhs -= rhs;
  return lhs;
}

/// @brief Subtracts a complex number from a real number
/// @param lhs Real number
/// @param rhs Complex number
/// @returns The difference of the two complex numbers
Complex operator-(const real_t lhs, const Complex& rhs)
{
  Complex tmp = -rhs;
  tmp += lhs;
  return tmp;
}

/// @brief Multiplies two complex numbers
/// @param lhs Complex number
/// @param rhs Complex number
/// @returns The multiplication of the two complex numbers
Complex operator*(Complex lhs, const Complex& rhs)
{
  lhs *= rhs;
  return lhs;
}

/// @brief Multiplies a real number and a complex number
/// @param lhs Complex number
/// @param rhs Real number
/// @returns The multiplication's result
Complex operator*(Complex lhs, const real_t rhs)
{
  lhs *= rhs;
  return lhs;
}

/// @brief Multiplies a complex number and a real number
/// @param lhs Real number
/// @param rhs Complex number
/// @returns The multiplication's result
Complex operator*(const real_t lhs, Complex rhs)
{
  rhs *= lhs;
  return rhs;
}

/// @brief Divides two complex numbers
/// @param lhs Complex number
/// @param rhs Complex number
/// @returns The division of the two complex numbers
Complex operator/(Complex lhs, const Complex& rhs)
{
  lhs /= rhs;
  return lhs;
}

/// @brief Divides a complex number by a real number
/// @param lhs Real number
/// @param rhs Complex number
/// @returns The division's result
Complex operator/(Complex lhs, const real_t rhs)
{
  lhs /= rhs;
  return lhs;
}

/// @brief Divides a real number by a complex number
/// @param lhs Real number
/// @param rhs Complex number
/// @returns The division's result
Complex operator/(const real_t lhs, const Complex& rhs)
{
  return Complex(lhs, 0) / rhs;
}

#ifdef __PRINT__

/// @brief Converts a complex number to a string
/// @param z Complex number
/// @returns The complex number as a string
std::string to_string(const Complex& z, const real_t epsilon)
{
  const std::string sign = z.imag() > 0 ? " + " : " - ";
  std::string str = std::to_string(z.real());
  if (std::abs(z.imag()) > epsilon) {
    str += sign + std::to_string(std::abs(z.imag())) + "i";
  }
  return str;
}

/// @brief Serializes a complex number
std::ostream& operator<<(std::ostream& os, const Complex& z)
{
  os << to_string(z);
  return os;
}

/// @brief Prints a complex number in the default output stream
void print(const Complex& z)
{
  std::cout << z << '\n';
}

#endif

}  // namespace mathlib

#ifdef __PRINT__
namespace std
{

/// @brief Overloads the conversion from a complex number to a string
std::string to_string(const mathlib::Complex& z)
{
  return mathlib::to_string(z);
}
#endif

}  // namespace std