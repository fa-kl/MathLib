/*****************************************************************************************
 * @file: Vector.hpp
 *
 * @brief: A class for n-dimensional vectors.
 *
 * @details: This file implements a class Vector for handling n-dimensional
 *vectors.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <memory>
#include <type_traits>
#include <vector>

#include "DivByZeroError.hpp"
#include "IncompatibleSizeError.hpp"
#include "IndexOutOfRangeError.hpp"
#include "MathLibError.hpp"
#include "compare.hpp"
#include "config.hpp"
#include "rounding.hpp"

#ifdef __PRINT__
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#endif

/*
  TODO:
    slice: [], ()
    sort
    unique
    rotate2d
    filter
    roll
    elem inversion
    dft
*/

namespace mathlib
{

#pragma region class Vector

template <typename data_t = real_t>
class Vector;

template <typename data_t>
class Matrix;

template <typename data_t>
class Vector
{
protected:
  std::unique_ptr<data_t[]> m_data;
  dim_t m_size;

public:
  Vector() : m_data(nullptr), m_size(0, 0) {}

  explicit Vector(size_t n) : m_data(std::make_unique<data_t[]>(n)), m_size(n, 1) {}

  Vector(std::initializer_list<data_t> values)
      : m_data(std::make_unique<data_t[]>(values.size())), m_size(values.size(), 1)
  {
    std::copy(values.begin(), values.end(), m_data.get());
  }

  Vector(const Vector& other)
      : m_data(other.length() > 0 ? std::make_unique<data_t[]>(other.length()) : nullptr), m_size(other.m_size)
  {
    if (m_data) {
      std::copy(other.m_data.get(), other.m_data.get() + length(), m_data.get());
    }
  }

  Vector(Vector&& other) noexcept : m_data(std::move(other.m_data)), m_size(other.m_size) { other.m_size = {0, 0}; }

  ~Vector() = default;

  Vector& operator=(const Vector& other)
  {
    if (this != &other) {
      if (other.length() > 0) {
        m_data = std::make_unique<data_t[]>(other.length());
        std::copy(other.m_data.get(), other.m_data.get() + other.length(), m_data.get());
      } else {
        m_data.reset();
      }
      m_size = other.m_size;
    }
    return *this;
  }

  Vector& operator=(Vector&& other) noexcept
  {
    if (this != &other) {
      m_data = std::move(other.m_data);
      m_size = other.m_size;
      other.m_size = {0, 0};
    }
    return *this;
  }

  size_t length() const { return m_size.rows * m_size.cols; }

  dim_t size() const { return m_size; }

  size_t rows() const { return m_size.rows; }

  size_t cols() const { return m_size.cols; }

  bool isEmpty() const { return length() == 0; }

  data_t* data() { return m_data.get(); }

  const data_t* data() const { return m_data.get(); }

  data_t& operator[](const size_t i)
  {
    if (i < length()) {
      return m_data[i];
    }
    throw IndexOutOfRangeError(i);
  }

  const data_t& operator[](const size_t i) const
  {
    if (i < length()) {
      return m_data[i];
    }
    throw IndexOutOfRangeError(i);
  }

  data_t& operator()(const index_t i)
  {
    if (i == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    if (static_cast<size_t>(std::abs(i)) <= length()) {
      return i > 0 ? m_data[static_cast<size_t>(i - 1)]
                   : m_data[static_cast<size_t>(static_cast<index_t>(length()) + i)];
    }
    throw IndexOutOfRangeError(i);
  }

  const data_t& operator()(const index_t i) const
  {
    if (i == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    if (static_cast<size_t>(std::abs(i)) <= length()) {
      return i > 0 ? m_data[static_cast<size_t>(i - 1)]
                   : m_data[static_cast<size_t>(static_cast<index_t>(length()) + i)];
    }
    throw IndexOutOfRangeError(i);
  }

  Vector operator+() const { return Vector(*this); }

  Vector operator-() const
  {
    Vector tmp;
    for (size_t i = 0; i < this->length(); ++i) {
      tmp.m_data[i] = -m_data[i];
    }
    return tmp;
  }
};

#pragma endregion

#pragma region Vector Constructor Methods

template <typename data_t>
Vector<data_t> zeros(size_t n)
{
  return Vector<data_t>(n);
}

template <typename data_t>
Vector<data_t> ones(size_t n)
{
  Vector<data_t> result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = static_cast<data_t>(1);
  }
  return result;
}

#pragma endregion

#pragma region Arithmetic Operators

template <typename data_t>
Vector<data_t> operator+(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] + rhs[i];
  }
  return result;
}

template <typename T1, typename T2>
auto operator+(const Vector<T1>& lhs,
               const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = decltype(std::declval<T1>() + std::declval<T2>());
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) + static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator+(const Vector<data_t>& lhs, const scalar_t& rhs)
{
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] + static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator+(const scalar_t& lhs, const Vector<data_t>& rhs)
{
  size_t len = rhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) + rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<data_t> operator-(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}

template <typename T1, typename T2>
auto operator-(const Vector<T1>& lhs,
               const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() - std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = decltype(std::declval<T1>() - std::declval<T2>());
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) - static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator-(const Vector<data_t>& lhs, const scalar_t& rhs)
{
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] - static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator-(const scalar_t& lhs, const Vector<data_t>& rhs)
{
  size_t len = rhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) - rhs[i];
  }
  return result;
}

template <typename data_t>
real_t operator*(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  real_t result = 0.0;
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result += static_cast<real_t>(lhs[i]) * static_cast<real_t>(rhs[i]);
  }
  return result;
}

template <typename T1, typename T2>
real_t operator*(const Vector<T1>& lhs, const Vector<T2>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  real_t result = 0.0;
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result += static_cast<real_t>(lhs[i]) * static_cast<real_t>(rhs[i]);
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator*(const Vector<data_t>& lhs, const scalar_t& rhs)
{
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] * static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator*(const scalar_t& lhs, const Vector<data_t>& rhs)
{
  size_t len = rhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) * rhs[i];
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator/(const Vector<data_t>& lhs, const scalar_t& rhs)
{
  if (isFuzzyEqual(static_cast<real_t>(rhs), 0.0)) {
    throw DivByZeroError();
  }
  size_t len = lhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] / static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Vector<data_t> operator/(const scalar_t& lhs, const Vector<data_t>& rhs)
{
  size_t len = rhs.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    if (isFuzzyEqual(static_cast<real_t>(rhs[i]), 0.0)) {
      throw DivByZeroError();
    }
    result[i] = static_cast<data_t>(lhs) / rhs[i];
  }
  return result;
}

#pragma endregion

#pragma region Logical Operators

inline Vector<bool> operator&&(const Vector<bool>& lhs, const Vector<bool>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] && rhs[i];
  }
  return result;
}

inline Vector<bool> operator||(const Vector<bool>& lhs, const Vector<bool>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] || rhs[i];
  }
  return result;
}

inline Vector<bool> operator^(const Vector<bool>& lhs, const Vector<bool>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] ^ rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator==(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] == rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator!=(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] != rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator>(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] > rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator>=(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] >= rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator<(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] < rhs[i];
  }
  return result;
}

template <typename data_t>
Vector<bool> operator<=(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] <= rhs[i];
  }
  return result;
}

#pragma endregion

#pragma region Methods

template <typename data_t>
size_t length(const Vector<data_t>& vec)
{
  return vec.length();
}

template <typename data_t>
dim_t size(const Vector<data_t>& vec)
{
  return vec.size();
}

template <typename data_t>
bool isEmpty(const Vector<data_t>& vec)
{
  return vec.isEmpty();
}

template <typename data_t>
Matrix<data_t> transpose(const Vector<data_t>& vec)
{
  Matrix<data_t> result(vec.cols(), vec.rows());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = vec[i];
  }
  return result;
}

/*
  TODO: cat, hcat, vcat, reshape, vec[i] extension, vectorize
  maybe also everything in Matrix!?
*/

template <typename data_t>
real_t norm(const Vector<data_t>& vec, real_t p = real_t(2))
{
  // +/- infinity norms
  if (std::isinf(p) && p > 0) {
    return max(abs(vec));
  }
  if (std::isinf(p) && p < 0) {
    return min(abs(vec));
  }
  // check for valid p
  if (p <= 0) {
    throw MathLibError("Invalid p");
  }
  // 1-norm
  if (p == 1.0) {
    return sum(abs(vec));
  }
  // 2-norm
  if (p == 2.0) {
    return std::sqrt(sum(pow(abs(vec), real_t(2))));
  }
  // p-norm
  return std::pow(sum(pow(abs(vec), p)), 1.0 / p);
}

template <typename data_t>
Vector<data_t> normalize(const Vector<data_t>& vec)
{
  real_t n = norm(vec);
  if (isFuzzyEqual(n, 0.0)) {
    throw MathLibError("Cannot normalize a zero-length vector.");
  }
  return vec / n;
}

template <typename data_t>
Vector<data_t> abs(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::abs(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> abs(const Vector<data_t>& vec);

template <typename data_t>
Vector<data_t> exp(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::exp(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> log(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::log(vec[i]);
  }
  return result;
}

#pragma region Trigonometric Methods

template <typename data_t>
Vector<data_t> sin(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::sin(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> cos(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::cos(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> tan(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::tan(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> asin(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::asin(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> acos(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::acos(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> atan(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::atan(vec[i]);
  }
  return result;
}

#pragma endregion

template <typename data_t>
Vector<data_t> pow(const Vector<data_t>& vec, real_t exponent)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::pow(vec[i], exponent);
  }
  return result;
}

template <typename data_t>
Vector<data_t> sqrt(const Vector<data_t>& vec)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    if (vec[i] < static_cast<data_t>(0)) {
      throw MathLibError("Square root of negative numbers not supported");
    }
    result[i] = std::sqrt(vec[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> round(const Vector<data_t>& vec, const size_t& decimal, RoundingMethod_t method = NEAREST)
{
  size_t len = vec.length();
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = round(vec[i], decimal, method);
  }
  return result;
}

template <typename data_t>
data_t sum(const Vector<data_t>& vec)
{
  data_t result = 0.0;
  size_t len = vec.length();
  for (size_t i = 0; i < len; ++i) {
    result += vec[i];
  }
  return result;
}

template <typename T1, typename T2>
auto sum(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = std::common_type_t<T1, T2>;
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) + static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t>
Vector<data_t> diff(const Vector<data_t>& vec)
{
  size_t len = vec.length() - 1;
  Vector<data_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = vec[i] - vec[i + 1];
  }
  return result;
}

template <typename T1, typename T2>
auto diff(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = std::common_type_t<T1, T2>;
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) + static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename T1, typename T2>
auto elmul(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = std::common_type_t<T1, T2>;
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) * static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename T1, typename T2>
real_t dot(const Vector<T1>& lhs, const Vector<T2>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  real_t result = 0.0;
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result += static_cast<real_t>(lhs[i]) * static_cast<real_t>(rhs[i]);
  }
  return result;
}

template <typename T1, typename T2>
auto cross(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != 3 || rhs.length() != 3) {
    throw IncompatibleSizeError("The cross product is only defined for 3-dimensional vectors");
  }
  using result_t = std::common_type_t<T1, T2>;
  return Vector<result_t>({static_cast<result_t>(lhs[1]) * static_cast<result_t>(rhs[2]) -
                               static_cast<result_t>(lhs[2]) * static_cast<result_t>(rhs[1]),
                           static_cast<result_t>(lhs[2]) * static_cast<result_t>(rhs[0]) -
                               static_cast<result_t>(lhs[0]) * static_cast<result_t>(rhs[2]),
                           static_cast<result_t>(lhs[0]) * static_cast<result_t>(rhs[1]) -
                               static_cast<result_t>(lhs[1]) * static_cast<result_t>(rhs[0])});
}

/*
  TODO: outer prod, hadamard_product
*/

template <typename data_t>
data_t max(const Vector<data_t>& vec)
{
  data_t result = -INF;
  size_t len = vec.length();
  for (size_t i = 0; i < len; ++i) {
    result = result > vec[i] ? result : vec[i];
  }
  return result;
}

template <typename data_t>
data_t min(const Vector<data_t>& vec)
{
  data_t result = INF;
  size_t len = vec.length();
  for (size_t i = 0; i < len; ++i) {
    result = result < vec[i] ? result : vec[i];
  }
  return result;
}

template <typename data_t>
real_t mean(const Vector<data_t>& vec)
{
  return sum(vec) / static_cast<real_t>(vec.length());
}

template <typename data_t>
real_t std(const Vector<data_t>& vec)
{
  return std::sqrt(var(vec));
}

template <typename data_t>
real_t var(const Vector<data_t>& vec)
{
  real_t factor = static_cast<real_t>(1.0) / (static_cast<real_t>(vec.length()) - static_cast<real_t>(1.0));
  Vector<real_t> v = abs(vec - mean(vec));
  return factor * sum(elmul(v, v));
}

// TODO: cov(const Vector<data_t>& vec)
// TODO: corr(const Vector<data_t>& vec)

inline bool all(const Vector<bool>& vec)
{
  size_t len = vec.length();
  for (size_t i = 0; i < len; ++i) {
    if (!vec[i]) {
      return false;
    }
  }
  return true;
}

inline bool any(const Vector<bool>& vec)
{
  size_t len = vec.length();
  for (size_t i = 0; i < len; ++i) {
    if (vec[i]) {
      return true;
    }
  }
  return false;
}

inline Vector<bool> isFuzzyEqual(const Vector<real_t>& lhs, const Vector<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzyEqual(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Vector<bool> isStrictFuzzyGreater(const Vector<real_t>& lhs, const Vector<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = isStrictFuzzyGreater(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Vector<bool> isFuzzyGreater(const Vector<real_t>& lhs, const Vector<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzyGreater(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Vector<bool> isStrictFuzzySmaller(const Vector<real_t>& lhs, const Vector<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = isStrictFuzzySmaller(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Vector<bool> isFuzzySmaller(const Vector<real_t>& lhs, const Vector<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  size_t len = lhs.length();
  Vector<bool> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzySmaller(lhs[i], rhs[i], epsilon);
  }
  return result;
}

#pragma endregion

#ifdef __PRINT__
template <typename data_t>
std::string to_string(const Vector<data_t>& vec, const size_t& decimal = 3, RoundingMethod_t method = NEAREST)
{
  size_t len = vec.length();
  if (len == 0) {
    return "[]";
  }

  // First pass: calculate maximum width needed
  size_t max_width = 0;
  for (size_t i = 0; i < len; ++i) {
    data_t rounded = mathlib::round(vec[i], decimal, method);
    std::ostringstream temp;
    temp << std::fixed << std::setprecision(static_cast<int>(decimal)) << rounded;
    max_width = std::max(max_width, temp.str().length());
  }

  // Second pass: format with alignment
  dim_t dim = vec.size();
  std::string seperator = dim.rows >= dim.cols ? "\n " : ", ";
  std::ostringstream oss;
  oss << "[";

  for (size_t i = 0; i < len; ++i) {
    data_t rounded = mathlib::round(vec[i], decimal, method);
    std::ostringstream temp;
    temp << std::fixed << std::setprecision(static_cast<int>(decimal)) << rounded;
    std::string value_str = temp.str();

    // Right-align the value
    oss << std::setw(static_cast<int>(max_width)) << value_str;

    if (i < len - 1) {
      oss << seperator;
    }
  }

  oss << "]";
  return oss.str();
}

template <typename data_t>
std::ostream& operator<<(std::ostream& os, const Vector<data_t>& vec)
{
  os << to_string(vec);
  return os;
}

template <typename data_t>
void print(const Vector<data_t>& vec)
{
  std::cout << vec << '\n';
}
#endif

}  // namespace mathlib

#ifdef __PRINT__

namespace std
{
template <typename data_t>
std::string to_string(const mathlib::Vector<data_t>& vec)
{
  return mathlib::to_string(vec);
}

#endif

}  // namespace std