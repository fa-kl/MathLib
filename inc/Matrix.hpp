/*****************************************************************************************
 * @file: Matrix.hpp
 *
 * @brief: A class for matrices.
 *
 * @details: This file implements a class Matrix for handling matrices.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>

#include "DivByZeroError.hpp"
#include "IncompatibleSizeError.hpp"
#include "IndexOutOfRangeError.hpp"
#include "MathLibError.hpp"
#include "Vector.hpp"
#include "compare.hpp"
#include "config.hpp"
#include "rounding.hpp"

#ifdef __PRINT__
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#endif

namespace mathlib
{

#pragma region class Matrix

template <typename data_t = real_t>
class Matrix;

template <typename data_t>
class Matrix
{
protected:
  std::unique_ptr<data_t[]> m_data;
  dim_t m_size;

public:
  Matrix() : m_data(nullptr), m_size({0, 0}) {}

  explicit Matrix(size_t rows, size_t cols) : m_data(std::make_unique<data_t[]>(rows * cols)), m_size({rows, cols}) {}

  Matrix(size_t rows, size_t cols, data_t value) : m_data(std::make_unique<data_t[]>(rows * cols)), m_size({rows, cols})
  {
    size_t len = length();
    for (size_t i = 0; i < len; ++i) {
      m_data[i] = value;
    }
  }

  Matrix(std::initializer_list<std::initializer_list<data_t>> values)
  {
    size_t rows = values.size();
    if (rows == 0) {
      m_data = nullptr;
      m_size = {0, 0};
      return;
    }
    size_t cols = values.begin()->size();
    for (const auto& row : values) {
      if (row.size() != cols) {
        throw MathLibError("All rows must have the same number of columns");
      }
    }
    m_data = std::make_unique<data_t[]>(rows * cols);
    m_size = {rows, cols};
    size_t row_idx = 0;
    for (const auto& row : values) {
      size_t col_idx = 0;
      for (const auto& element : row) {
        m_data[col_idx * rows + row_idx] = element;
        ++col_idx;
      }
      ++row_idx;
    }
  }

  Matrix(std::initializer_list<Vector<data_t>> vectors)
  {
    size_t cols = vectors.size();
    if (cols == 0) {
      m_data = nullptr;
      m_size = {0, 0};
      return;
    }
    size_t rows = vectors.begin()->length();
    for (const auto& vec : vectors) {
      if (vec.length() != rows) {
        throw MathLibError("All vectors must have the same length");
      }
    }
    m_data = std::make_unique<data_t[]>(rows * cols);
    m_size = {rows, cols};
    size_t col_idx = 0;
    for (const auto& vec : vectors) {
      for (size_t row_idx = 0; row_idx < rows; ++row_idx) {
        m_data[col_idx * rows + row_idx] = vec[row_idx];
      }
      ++col_idx;
    }
  }

  explicit Matrix(const Vector<data_t>& vec) : m_data(std::make_unique<data_t[]>(vec.length())), m_size(vec.size())
  {
    for (size_t i = 0; i < vec.length(); ++i) {
      m_data[i] = vec[i];
    }
  }

  Matrix(const Matrix& other)
      : m_data(other.length() > 0 ? std::make_unique<data_t[]>(other.length()) : nullptr), m_size(other.m_size)
  {
    if (m_data) {
      std::copy(other.m_data.get(), other.m_data.get() + length(), m_data.get());
    }
  }

  Matrix(Matrix&& other) noexcept : m_data(std::move(other.m_data)), m_size(other.m_size) { other.m_size = {0, 0}; }

  ~Matrix() = default;

  Matrix& operator=(const Matrix& other)
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

  Matrix& operator=(Matrix&& other) noexcept
  {
    if (this != &other) {
      m_data = std::move(other.m_data);
      m_size = other.m_size;
      other.m_size = {0, 0};
    }
    return *this;
  }

  dim_t size() const { return m_size; }

  size_t length() const { return m_size.rows * m_size.cols; }

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
    throw IndexOutOfRangeError();
  }

  const data_t& operator[](const size_t i) const
  {
    if (i < length()) {
      return m_data[i];
    }
    throw IndexOutOfRangeError();
  }

  template <typename int_t, typename = std::enable_if_t<std::is_integral_v<int_t>>>
  data_t& operator()(const int_t i)
  {
    if (i == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    index_t idx_val = static_cast<index_t>(i);
    index_t len = static_cast<index_t>(length());
    size_t idx = idx_val > 0 ? static_cast<size_t>(idx_val - 1) : static_cast<size_t>(len + idx_val);
    if (idx >= length()) {
      throw IndexOutOfRangeError();
    }
    return m_data[idx];
  }

  template <typename int_t, typename = std::enable_if_t<std::is_integral_v<int_t>>>
  const data_t& operator()(const int_t i) const
  {
    if (i == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    index_t idx_val = static_cast<index_t>(i);
    index_t len = static_cast<index_t>(length());
    size_t idx = idx_val > 0 ? static_cast<size_t>(idx_val - 1) : static_cast<size_t>(len + idx_val);
    if (idx >= length()) {
      throw IndexOutOfRangeError();
    }
    return m_data[idx];
  }

  template <typename int_t1,
            typename int_t2,
            typename = std::enable_if_t<std::is_integral_v<int_t1> && std::is_integral_v<int_t2>>>
  data_t& operator()(const int_t1 row, const int_t2 col)
  {
    if (row == 0 || col == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    index_t row_val = static_cast<index_t>(row);
    index_t col_val = static_cast<index_t>(col);
    size_t r = row_val > 0 ? static_cast<size_t>(row_val - 1)
                           : static_cast<size_t>(static_cast<index_t>(m_size.rows) + row_val);
    size_t c = col_val > 0 ? static_cast<size_t>(col_val - 1)
                           : static_cast<size_t>(static_cast<index_t>(m_size.cols) + col_val);
    if (r >= m_size.rows || c >= m_size.cols) {
      throw IndexOutOfRangeError();
    }
    return m_data[c * m_size.rows + r];
  }

  template <typename int_t1,
            typename int_t2,
            typename = std::enable_if_t<std::is_integral_v<int_t1> && std::is_integral_v<int_t2>>>
  const data_t& operator()(const int_t1 row, const int_t2 col) const
  {
    if (row == 0 || col == 0) {
      throw IndexOutOfRangeError("0 is an invalid index for 1-based access methods");
    }
    index_t row_val = static_cast<index_t>(row);
    index_t col_val = static_cast<index_t>(col);
    size_t r = row_val > 0 ? static_cast<size_t>(row_val - 1)
                           : static_cast<size_t>(static_cast<index_t>(m_size.rows) + row_val);
    size_t c = col_val > 0 ? static_cast<size_t>(col_val - 1)
                           : static_cast<size_t>(static_cast<index_t>(m_size.cols) + col_val);
    if (r >= m_size.rows || c >= m_size.cols) {
      throw IndexOutOfRangeError();
    }
    return m_data[c * m_size.rows + r];
  }

  Matrix operator+() const { return Matrix(*this); }

  Matrix operator-() const
  {
    Matrix result(m_size.rows, m_size.cols);
    size_t len = length();
    for (size_t i = 0; i < len; ++i) {
      result.m_data[i] = -m_data[i];
    }
    return result;
  }
};

#pragma endregion

#pragma region Matrix Constructor Methods

template <typename data_t = real_t>
Matrix<data_t> zeros(size_t rows, size_t cols)
{
  return Matrix<data_t>(rows, cols);
}

template <typename data_t = real_t>
Matrix<data_t> ones(size_t rows, size_t cols)
{
  return Matrix<data_t>(rows, cols, static_cast<data_t>(1));
}

template <typename data_t = real_t>
Matrix<data_t> eye(size_t n)
{
  Matrix<data_t> result(n, n);
  for (size_t i = 0; i < n; ++i) {
    result[i * n + i] = static_cast<data_t>(1);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> diag(std::initializer_list<data_t> values)
{
  size_t n = values.size();
  Matrix<data_t> result(n, n);
  index_t i = 1;
  for (const auto& val : values) {
    result(i, i) = val;
    i++;
  }
  return result;
}

template <typename data_t>
Matrix<data_t> diag(const Vector<data_t>& vec)
{
  size_t n = vec.length();
  Matrix<data_t> result(n, n);
  for (size_t i = 0; i < n; ++i) {
    result[i * n + i] = vec[i];
  }
  return result;
}

template <typename data_t>
Vector<data_t> diag(const Matrix<data_t>& mat)
{
  size_t n = std::min(mat.rows(), mat.cols());
  Vector<data_t> result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = mat[i * mat.rows() + i];
  }
  return result;
}

#pragma endregion

#pragma region Arithmetic Operators

template <typename data_t>
Matrix<data_t> operator+(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] + rhs[i];
  }
  return result;
}

template <typename T1, typename T2>
auto operator+(const Matrix<T1>& lhs,
               const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  using result_t = decltype(std::declval<T1>() + std::declval<T2>());
  Matrix<result_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) + static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator+(const Matrix<data_t>& lhs, const scalar_t& rhs)
{
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] + static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator+(const scalar_t& lhs, const Matrix<data_t>& rhs)
{
  Matrix<data_t> result(rhs.rows(), rhs.cols());
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) + rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> operator-(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}

template <typename T1, typename T2>
auto operator-(const Matrix<T1>& lhs,
               const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() - std::declval<T2>())>
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  using result_t = decltype(std::declval<T1>() - std::declval<T2>());
  Matrix<result_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) - static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator-(const Matrix<data_t>& lhs, const scalar_t& rhs)
{
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] - static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator-(const scalar_t& lhs, const Matrix<data_t>& rhs)
{
  Matrix<data_t> result(rhs.rows(), rhs.cols());
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) - rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> operator*(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.cols() != rhs.rows()) {
    throw IncompatibleSizeError();
  }
  Matrix<data_t> result(lhs.rows(), rhs.cols());
  for (size_t i = 0; i < lhs.rows(); ++i) {
    for (size_t col = 0; col < rhs.cols(); ++col) {
      data_t sum = static_cast<data_t>(0);
      for (size_t k = 0; k < lhs.cols(); ++k) {
        sum += lhs[k * lhs.rows() + i] * rhs[col * rhs.rows() + k];
      }
      result[col * result.rows() + i] = sum;
    }
  }
  return result;
}

template <typename T1, typename T2>
auto operator*(const Matrix<T1>& lhs,
               const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() * std::declval<T2>())>
{
  if (lhs.cols() != rhs.rows()) {
    throw IncompatibleSizeError();
  }
  using result_t = decltype(std::declval<T1>() * std::declval<T2>());
  Matrix<result_t> result(lhs.rows(), rhs.cols());
  for (size_t i = 0; i < lhs.rows(); ++i) {
    for (size_t col = 0; col < rhs.cols(); ++col) {
      result_t sum = static_cast<result_t>(0);
      for (size_t k = 0; k < lhs.cols(); ++k) {
        sum += static_cast<result_t>(lhs[k * lhs.rows() + i]) * static_cast<result_t>(rhs[col * rhs.rows() + k]);
      }
      result[col * result.rows() + i] = sum;
    }
  }
  return result;
}

template <typename data_t>
Vector<data_t> operator*(const Matrix<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.cols() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  Vector<data_t> result(lhs.rows());
  for (size_t i = 0; i < lhs.rows(); ++i) {
    data_t sum = static_cast<data_t>(0);
    for (size_t col = 0; col < lhs.cols(); ++col) {
      sum += lhs[col * lhs.rows() + i] * rhs[col];
    }
    result[i] = sum;
  }
  return result;
}

template <typename T1, typename T2>
auto operator*(const Matrix<T1>& lhs,
               const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() * std::declval<T2>())>
{
  if (lhs.cols() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = decltype(std::declval<T1>() * std::declval<T2>());
  Vector<result_t> result(lhs.rows());
  for (size_t i = 0; i < lhs.rows(); ++i) {
    result_t sum = static_cast<result_t>(0);
    for (size_t col = 0; col < lhs.cols(); ++col) {
      sum += static_cast<result_t>(lhs[col * lhs.rows() + i]) * static_cast<result_t>(rhs[col]);
    }
    result[i] = sum;
  }
  return result;
}

template <typename data_t>
Matrix<data_t> operator*(const Vector<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (1 != rhs.rows()) {
    throw IncompatibleSizeError();
  }
  Matrix<data_t> result(lhs.length(), rhs.cols());
  for (size_t i = 0; i < lhs.length(); ++i) {
    for (size_t col = 0; col < rhs.cols(); ++col) {
      result[col * result.rows() + i] = lhs[i] * rhs[col];
    }
  }
  return result;
}

template <typename T1, typename T2>
auto operator*(const Vector<T1>& lhs,
               const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() * std::declval<T2>())>
{
  if (1 != rhs.rows()) {
    throw IncompatibleSizeError();
  }
  using result_t = decltype(std::declval<T1>() * std::declval<T2>());
  Matrix<result_t> result(lhs.length(), rhs.cols());
  for (size_t i = 0; i < lhs.length(); ++i) {
    for (size_t col = 0; col < rhs.cols(); ++col) {
      result[col * result.rows() + i] = static_cast<result_t>(lhs[i]) * static_cast<result_t>(rhs[col]);
    }
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator*(const Matrix<data_t>& lhs, const scalar_t& rhs)
{
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] * static_cast<data_t>(rhs);
  }
  return result;
}

template <typename scalar_t, typename data_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator*(const scalar_t& lhs, const Matrix<data_t>& rhs)
{
  Matrix<data_t> result(rhs.rows(), rhs.cols());
  size_t len = rhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<data_t>(lhs) * rhs[i];
  }
  return result;
}

template <typename data_t, typename scalar_t, typename = std::enable_if_t<is_scalar_type_v<scalar_t>>>
Matrix<data_t> operator/(const Matrix<data_t>& lhs, const scalar_t& rhs)
{
  if (isFuzzyEqual(static_cast<real_t>(rhs), 0.0)) {
    throw DivByZeroError();
  }
  Matrix<data_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] / static_cast<data_t>(rhs);
  }
  return result;
}

#pragma endregion

#pragma region Logical Operators

template <typename data_t>
Matrix<bool> operator==(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] == rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<bool> operator!=(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] != rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<bool> operator>(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] > rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<bool> operator>=(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] >= rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<bool> operator<(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] < rhs[i];
  }
  return result;
}

template <typename data_t>
Matrix<bool> operator<=(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] <= rhs[i];
  }
  return result;
}

inline Matrix<bool> operator&&(const Matrix<bool>& lhs, const Matrix<bool>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] && rhs[i];
  }
  return result;
}

inline Matrix<bool> operator||(const Matrix<bool>& lhs, const Matrix<bool>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] || rhs[i];
  }
  return result;
}

inline Matrix<bool> operator^(const Matrix<bool>& lhs, const Matrix<bool>& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = lhs[i] ^ rhs[i];
  }
  return result;
}

#pragma endregion

#pragma region Methods

template <typename data_t>
size_t rows(const Matrix<data_t>& mat)
{
  return mat.rows();
}

template <typename data_t>
size_t cols(const Matrix<data_t>& mat)
{
  return mat.cols();
}

template <typename data_t>
dim_t size(const Matrix<data_t>& mat)
{
  return mat.size();
}

template <typename data_t>
size_t length(const Matrix<data_t>& mat)
{
  return mat.length();
}

template <typename data_t>
bool isEmpty(const Matrix<data_t>& mat)
{
  return mat.isEmpty();
}

template <typename data_t>
Matrix<data_t> transpose(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.cols(), mat.rows());
  for (size_t i = 0; i < mat.rows(); ++i) {
    for (size_t col = 0; col < mat.cols(); ++col) {
      result[i * result.rows() + col] = mat[col * mat.rows() + i];
    }
  }
  return result;
}

template <typename data_t>
data_t trace(const Matrix<data_t>& mat)
{
  size_t n = std::min(mat.rows(), mat.cols());
  data_t result = static_cast<data_t>(0);
  for (size_t i = 0; i < n; ++i) {
    result += mat[i * mat.rows() + i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> abs(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::abs(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> exp(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::exp(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> log(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::log(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> sin(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::sin(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> cos(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::cos(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> tan(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::tan(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> asin(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::asin(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> acos(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::acos(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> atan(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::atan(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> elpow(const Matrix<data_t>& mat, real_t exponent)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::pow(mat[i], exponent);
  }
  return result;
}

// TODO:
/*
template <typename data_t>
Matrix<data_t> pow(const Matrix<data_t>& mat, real_t exponent)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = std::pow(mat[i], exponent);
  }
  return result;
}
*/

template <typename data_t>
Matrix<data_t> sqrt(const Matrix<data_t>& mat)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    if (mat[i] < static_cast<data_t>(0)) {
      throw MathLibError("Square root of negative numbers not supported");
    }
    result[i] = std::sqrt(mat[i]);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> round(const Matrix<data_t>& mat, const size_t& decimal, RoundingMethod_t method = NEAREST)
{
  Matrix<data_t> result(mat.rows(), mat.cols());
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = round(mat[i], decimal, method);
  }
  return result;
}

template <typename data_t>
data_t sum(const Matrix<data_t>& mat)
{
  data_t result = static_cast<data_t>(0);
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result += mat[i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> sum(const Matrix<data_t>& mat, const size_t dim)
{
  if (dim == 1) {
    Matrix<data_t> result(1, mat.cols());
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      data_t rowSum = 0.0;
      for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
        rowSum += mat(r, c);
      }
      result(c) = rowSum;
    }
    return result;
  }
  if (dim == 2) {
    Matrix<data_t> result(mat.rows(), 1);
    for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
      data_t colSum = 0.0;
      for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
        colSum += mat(r, c);
      }
      result(r) = colSum;
    }
    return result;
  }
  throw MathLibError("dim can only be 1 or 2");
}

template <typename T1, typename T2>
auto elmul(const Matrix<T1>& lhs, const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() * std::declval<T2>())>
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  using result_t = decltype(std::declval<T1>() * std::declval<T2>());
  Matrix<result_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) * static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename T1, typename T2>
auto eldiv(const Matrix<T1>& lhs, const Matrix<T2>& rhs) -> Matrix<decltype(std::declval<T1>() / std::declval<T2>())>
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  using result_t = decltype(std::declval<T1>() / std::declval<T2>());
  Matrix<result_t> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    if (isFuzzyEqual(static_cast<real_t>(rhs[i]), 0.0)) {
      throw DivByZeroError();
    }
    result[i] = static_cast<result_t>(lhs[i]) / static_cast<result_t>(rhs[i]);
  }
  return result;
}

template <typename data_t>
data_t max(const Matrix<data_t>& mat)
{
  data_t result = -INF;
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result = result > mat[i] ? result : mat[i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> max(const Matrix<data_t>& mat, const size_t dim)
{
  if (dim == 1) {
    Matrix<data_t> result(1, mat.cols());
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      data_t rowMax = -INF;
      for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
        rowMax = rowMax > mat(r, c) ? rowMax : mat(r, c);
      }
      result(c) = rowMax;
    }
    return result;
  }
  if (dim == 2) {
    Matrix<data_t> result(mat.rows(), 1);
    for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
      data_t colMax = -INF;
      for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
        colMax = colMax > mat(r, c) ? colMax : mat(r, c);
      }
      result(r) = colMax;
    }
    return result;
  }
  throw MathLibError("dim can only be 1 or 2");
}

template <typename data_t>
data_t min(const Matrix<data_t>& mat)
{
  data_t result = INF;
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    result = result < mat[i] ? result : mat[i];
  }
  return result;
}

template <typename data_t>
Matrix<data_t> min(const Matrix<data_t>& mat, const size_t dim)
{
  if (dim == 1) {
    Matrix<data_t> result(1, mat.cols());
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      data_t rowMin = INF;
      for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
        rowMin = rowMin < mat(r, c) ? rowMin : mat(r, c);
      }
      result(c) = rowMin;
    }
    return result;
  }
  if (dim == 2) {
    Matrix<data_t> result(mat.rows(), 1);
    for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
      data_t colMin = INF;
      for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
        colMin = colMin < mat(r, c) ? colMin : mat(r, c);
      }
      result(r) = colMin;
    }
    return result;
  }
  throw MathLibError("dim can only be 1 or 2");
}

template <typename data_t>
real_t mean(const Matrix<data_t>& mat)
{
  if (mat.length() == 0) {
    throw MathLibError("Cannot compute mean of empty matrix");
  }
  return static_cast<real_t>(sum(mat)) / static_cast<real_t>(mat.length());
}

template <typename data_t>
Matrix<real_t> mean(const Matrix<data_t>& mat, const size_t dim)
{
  if (mat.length() == 0) {
    throw MathLibError("Cannot compute mean of empty matrix");
  }
  if (dim == 1) {
    return sum(mat, dim) / static_cast<real_t>(mat.rows());
  }
  if (dim == 2) {
    return sum(mat, dim) / static_cast<real_t>(mat.cols());
  }
  throw MathLibError("dim can only be 1 or 2");
}

template <typename data_t>
real_t std(const Matrix<data_t>& mat)
{
  if (mat.length() == 0) {
    throw MathLibError("Cannot compute standard deviation of empty matrix");
  }
  return std::sqrt(var(mat));
}

template <typename data_t>
Matrix<real_t> std(const Matrix<data_t>& mat, const size_t dim)
{
  if (mat.length() == 0) {
    throw MathLibError("Cannot compute standard deviation of empty matrix");
  }
  return sqrt(var(mat, dim));
}

template <typename data_t>
real_t var(const Matrix<data_t>& mat)
{
  if (mat.length() == 0) {
    throw MathLibError("Cannot compute variance of empty matrix");
  }
  real_t m = mean(mat);
  real_t sum_sq_diff = 0.0;
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    real_t diff = static_cast<real_t>(mat[i]) - m;
    sum_sq_diff += diff * diff;
  }
  return sum_sq_diff / static_cast<real_t>(len);
}

template <typename data_t>
Matrix<real_t> var(const Matrix<data_t>& mat, const size_t dim)
{
  if (dim == 1) {
    Matrix<real_t> means = mean(mat, dim);
    Matrix<real_t> result(1, mat.cols());    
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      real_t sum_sq_diff = 0.0;
      for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
        real_t diff = static_cast<real_t>(mat(r, c)) - means(c);
        sum_sq_diff += diff * diff;
      }
      result(c) = sum_sq_diff / static_cast<real_t>(mat.rows());
    }
    return result;
  }
  if (dim == 2) {
    Matrix<real_t> means = mean(mat, dim);
    Matrix<real_t> result(mat.rows(), 1);    
    for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
      real_t sum_sq_diff = 0.0;
      for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
        real_t diff = static_cast<real_t>(mat(r, c)) - means(r);
        sum_sq_diff += diff * diff;
      }
      result(r) = sum_sq_diff / static_cast<real_t>(mat.cols());
    }
    return result;
  }
  throw MathLibError("dim can only be 1 or 2");
}

inline bool all(const Matrix<bool>& mat)
{
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    if (!mat[i]) {
      return false;
    }
  }
  return true;
}

inline bool any(const Matrix<bool>& mat)
{
  size_t len = mat.length();
  for (size_t i = 0; i < len; ++i) {
    if (mat[i]) {
      return true;
    }
  }
  return false;
}

template <typename data_t>
Matrix<data_t> hcat(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.rows() != rhs.rows()) {
    throw IncompatibleSizeError("Matrices must have the same number of rows for horizontal concatenation");
  }
  Matrix<data_t> result(lhs.rows(), lhs.cols() + rhs.cols());
  for (size_t col = 1; col <= lhs.cols(); ++col) {
    for (size_t row = 1; row <= lhs.rows(); ++row) {
      result(row, col) = lhs(row, col);
    }
  }
  for (size_t col = 1; col <= rhs.cols(); ++col) {
    for (size_t row = 1; row <= rhs.rows(); ++row) {
      result(row, lhs.cols() + col) = rhs(row, col);
    }
  }
  return result;
}

template <typename data_t>
Matrix<data_t> vcat(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError("Matrices must have the same number of columns for vertical concatenation");
  }
  Matrix<data_t> result(lhs.rows() + rhs.rows(), lhs.cols());
  for (index_t col = 1; col <= static_cast<index_t>(lhs.cols()); ++col) {
    for (index_t row = 1; row <= static_cast<index_t>(lhs.rows()); ++row) {
      result(row, col) = lhs(row, col);
    }
  }
  for (index_t col = 1; col <= static_cast<index_t>(rhs.cols()); ++col) {
    for (index_t row = 1; row <= static_cast<index_t>(rhs.rows()); ++row) {
      result(static_cast<index_t>(lhs.rows()) + row, col) = rhs(row, col);
    }
  }
  return result;
}

template <typename data_t>
Matrix<data_t> hcat(const Matrix<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.rows() != rhs.length()) {
    throw IncompatibleSizeError("Matrix and vector must have compatible dimensions for horizontal concatenation");
  }
  Matrix<data_t> result(lhs.rows(), lhs.cols() + 1);
  for (index_t col = 1; col <= static_cast<index_t>(lhs.cols()); ++col) {
    for (index_t row = 1; row <= static_cast<index_t>(lhs.rows()); ++row) {
      result(row, col) = lhs(row, col);
    }
  }
  index_t last_col = static_cast<index_t>(lhs.cols()) + 1;
  for (index_t row = 1; row <= static_cast<index_t>(rhs.length()); ++row) {
    result(row, last_col) = rhs(row);
  }

  return result;
}

template <typename data_t>
Matrix<data_t> hcat(const Vector<data_t>& lhs, const Matrix<data_t>& rhs)
{
  if (lhs.length() != rhs.rows()) {
    throw IncompatibleSizeError("Vector and matrix must have compatible dimensions for horizontal concatenation");
  }
  Matrix<data_t> result(lhs.length(), 1 + rhs.cols());
  for (index_t row = 1; row <= static_cast<index_t>(lhs.length()); ++row) {
    result(row, 1) = lhs(row);
  }
  for (index_t col = 1; col <= static_cast<index_t>(rhs.cols()); ++col) {
    for (index_t row = 1; row <= static_cast<index_t>(rhs.rows()); ++row) {
      result(row, col + 1) = rhs(row, col);
    }
  }

  return result;
}

template <typename data_t>
Matrix<data_t> cat(const Matrix<data_t>& lhs, const Matrix<data_t>& rhs, size_t dim)
{
  if (dim == 1) {
    return vcat(lhs, rhs);
  } else if (dim == 2) {
    return hcat(lhs, rhs);
  } else {
    throw MathLibError("cat only supports dim=1 (vertical) or dim=2 (horizontal)");
  }
}

template <typename data_t>
real_t norm(const Matrix<data_t>& mat, real_t p = 2.0)
{
  return std::sqrt(static_cast<real_t>(trace(transpose(mat) * mat)));
}

template <typename data_t>
Matrix<data_t> vecnorm(const Matrix<data_t>& mat, real_t p = 2.0, size_t dim = 1)
{
  // +/- infinity norms
  if (std::isinf(p) && p > 0) {
    return max(abs(mat), dim);
  }
  if (std::isinf(p) && p < 0) {
    return min(abs(mat), dim);
  }
  // check for valid p
  if (p <= 0) {
    throw MathLibError("Invalid p");
  }
  // 1-norm
  if (p == 1.0) {
    return sum(abs(mat), dim);
  }
  // 2-norm
  if (p == 2.0) {
    return sqrt(sum(elpow(abs(mat), p), dim));
  }
  // p-norm
  return elpow(sum(elpow(abs(mat), p), dim), 1.0 / p);
}

template <typename data_t>
Matrix<real_t> dot(const Matrix<data_t>& A, const Matrix<data_t>& B, size_t dim = 1)
{
  if (A.rows() != B.rows() || A.cols() != B.cols()) {
    throw IncompatibleSizeError(A.rows(), A.cols(), B.rows(), B.cols());
  }

  if (dim == 1) {
    // Dot product along columns (result has one element per column)
    size_t n = A.cols();
    Matrix<real_t> result(1, n);

    for (size_t col = 0; col < n; ++col) {
      real_t sum = 0.0;
      for (size_t i = 0; i < A.rows(); ++i) {
        sum += static_cast<real_t>(A[col * A.rows() + i]) * static_cast<real_t>(B[col * B.rows() + i]);
      }
      result[col] = sum;
    }
    return result;
  } else if (dim == 2) {
    // Dot product along rows (result has one element per row)
    size_t m = A.rows();
    Matrix<real_t> result(m, 1);

    for (size_t row = 0; row < m; ++row) {
      real_t sum = 0.0;
      for (size_t col = 0; col < A.cols(); ++col) {
        sum += static_cast<real_t>(A[col * A.rows() + row]) * static_cast<real_t>(B[col * B.rows() + row]);
      }
      result[row] = sum;
    }
    return result;
  } else {
    throw MathLibError("dot only supports dim=1 (columns) or dim=2 (rows)");
  }
}

template <typename data_t>
Matrix<data_t> normalize(const Matrix<data_t>& mat, int p = 2, size_t dim = 1)
{
  if (dim == 1) {
    // Normalize each column
    Matrix<data_t> result(mat);
    Matrix<real_t> norms = vecnorm(mat, p, 1);

    for (size_t col = 0; col < mat.cols(); ++col) {
      if (isFuzzyEqual(norms[col], 0.0)) {
        throw MathLibError("Cannot normalize zero column");
      }
      for (size_t i = 0; i < mat.rows(); ++i) {
        result[col * mat.rows() + i] /= static_cast<data_t>(norms[col]);
      }
    }
    return result;
  } else if (dim == 2) {
    // Normalize each row
    Matrix<data_t> result(mat);
    Matrix<real_t> norms = vecnorm(mat, p, 2);

    for (size_t row = 0; row < mat.rows(); ++row) {
      if (isFuzzyEqual(norms[row], 0.0)) {
        throw MathLibError("Cannot normalize zero row");
      }
      for (size_t col = 0; col < mat.cols(); ++col) {
        result[col * mat.rows() + row] /= static_cast<data_t>(norms[row]);
      }
    }
    return result;
  } else {
    throw MathLibError("normalize only supports dim=1 (columns) or dim=2 (rows)");
  }
}

template <typename data_t>
data_t colDot(const Matrix<data_t>& mat, index_t col1, index_t col2)
{
  if (static_cast<size_t>(std::abs(col1)) > mat.cols()) {
    throw IndexOutOfRangeError(col1);
  }
  if (static_cast<size_t>(std::abs(col2)) > mat.cols()) {
    throw IndexOutOfRangeError(col2);
  }
  data_t result = 0;
  for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
    result += mat(r, col1) * mat(r, col2);
  }
  return result;
}

template <typename data_t>
data_t rowDot(const Matrix<data_t>& mat, index_t row1, index_t row2)
{
  if (static_cast<size_t>(std::abs(row1)) > mat.rows()) {
    throw IndexOutOfRangeError(row1);
  }
  if (static_cast<size_t>(std::abs(row2)) > mat.rows()) {
    throw IndexOutOfRangeError(row2);
  }
  data_t result = 0;
  for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
    result += mat(row1, c) * mat(row2, c);
  }
  return result;
}

template <typename data_t>
Matrix<data_t> orth(const Matrix<data_t>& A, real_t tol = EPSILON)
{
  size_t m = A.rows();
  size_t n = A.cols();

  if (m < n) {
    throw MathLibError("Orthogonalization requires rows >= cols");
  }

  // Copy matrix for orthogonalization
  Matrix<data_t> Q(A);
  std::vector<index_t> valid_cols;
  valid_cols.reserve(n);

  // Modified Gram-Schmidt Process
  // Optimization: Uses column-oriented operations for better cache locality in column-major storage
  for (index_t col = 1; col <= static_cast<index_t>(n); ++col) {
    // Orthogonalize against all previous orthonormal columns
    for (index_t prev : valid_cols) {
      // Compute dot product: <q_prev, q_col>
      data_t dot = static_cast<data_t>(0);
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        dot += Q(row, prev) * Q(row, col);
      }

      // Subtract projection: q_col = q_col - dot * q_prev
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        Q(row, col) -= dot * Q(row, prev);
      }
    }

    // Compute norm of remaining component
    // Note: vecnorm returns Matrix, need to extract single element
    data_t norm_squared = static_cast<data_t>(0);
    for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
      norm_squared += Q(row, col) * Q(row, col);
    }
    real_t norm = std::sqrt(static_cast<real_t>(norm_squared));

    // Only keep column if it's not zero (within tolerance)
    if (norm > tol) {
      // Normalize the column
      data_t norm_inv = static_cast<data_t>(1.0) / static_cast<data_t>(norm);
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        Q(row, col) *= norm_inv;
      }
      valid_cols.push_back(col);
    }
  }

  // Build result matrix with only valid orthonormal columns
  size_t rank = valid_cols.size();
  Matrix<data_t> result(m, rank);

  for (index_t new_col = 1; new_col <= static_cast<index_t>(rank); ++new_col) {
    index_t old_col = valid_cols[static_cast<size_t>(new_col - 1)];
    for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
      result(row, new_col) = Q(row, old_col);
    }
  }

  return result;
}

template <typename data_t>
size_t rank(const Matrix<data_t>& mat, real_t tol = EPSILON)
{
  if (mat.isEmpty()) {
    return 0;
  }

  size_t m = mat.rows();
  size_t n = mat.cols();

  // Use QR decomposition approach: rank = number of non-zero diagonal elements in R
  // Handle the case where m < n by transposing
  size_t min_dim = std::min(m, n);

  Matrix<data_t> A = m >= n ? mat : transpose(mat);
  Matrix<data_t> Q(A);
  size_t r = 0;

  // Modified Gram-Schmidt to count linearly independent columns
  // Optimization: Process only up to min_dim, use multiplication instead of division
  for (index_t col = 1; col <= static_cast<index_t>(min_dim); ++col) {
    // Orthogonalize against all previous orthonormal columns
    for (index_t prev = 1; prev <= static_cast<index_t>(r); ++prev) {
      // Compute dot product: <q_prev, q_col>
      data_t dot = static_cast<data_t>(0);
      for (index_t row = 1; row <= static_cast<index_t>(A.rows()); ++row) {
        dot += Q(row, prev) * Q(row, col);
      }

      // Subtract projection: q_col = q_col - dot * q_prev
      for (index_t row = 1; row <= static_cast<index_t>(A.rows()); ++row) {
        Q(row, col) -= dot * Q(row, prev);
      }
    }

    // Check if column is non-zero (linearly independent)
    // Compute norm squared directly for efficiency
    data_t norm_squared = static_cast<data_t>(0);
    for (index_t row = 1; row <= static_cast<index_t>(A.rows()); ++row) {
      norm_squared += Q(row, col) * Q(row, col);
    }
    real_t norm = std::sqrt(static_cast<real_t>(norm_squared));

    if (norm > tol) {
      // Normalize for next iteration (multiply by inverse to avoid repeated divisions)
      data_t norm_inv = static_cast<data_t>(1.0) / static_cast<data_t>(norm);
      for (index_t row = 1; row <= static_cast<index_t>(A.rows()); ++row) {
        Q(row, col) *= norm_inv;
      }
      r++;
    }
  }

  return r;
}

template <typename data_t>
bool isFullRank(const Matrix<data_t>& mat, real_t tol = EPSILON)
{
  return rank(mat, tol) == std::min(mat.rows(), mat.cols());
}

struct QRDecomposition {
  Matrix<real_t> Q;
  Matrix<real_t> R;
};

template<typename data_t>
QRDecomposition qr(const Matrix<data_t>& A, bool checkRank = true)
{
  size_t m = A.rows();
  size_t n = A.cols();
  if (m < n) {
    throw MathLibError("QR decomposition requires rows >= cols");
  }
  // Initialize Q as copy of A (will be transformed in place)
  // Initialize R as zeros
  Matrix<real_t> Q(A);
  Matrix<real_t> R(n, n);
  // Modified Gram-Schmidt Process
  // Optimization: Column-oriented operations for cache-friendly access pattern
  // Optimization: Use multiplication by inverse instead of repeated division
  for (index_t col = 1; col <= static_cast<index_t>(n); ++col) {
    // Orthogonalize column col against all previous orthonormal columns
    for (index_t prev = 1; prev < col; ++prev) {
      // R[prev,col] = <q_prev, a_col> = dot product
      real_t dot = static_cast<real_t>(0);
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        dot += Q(row, prev) * Q(row, col);
      }
      // Store in R matrix
      R(prev, col) = dot;
      // Remove component: a_col = a_col - R[prev,col] * q_prev
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        Q(row, col) -= dot * Q(row, prev);
      }
    }
    // Compute norm of remaining component
    // Direct computation avoids function call overhead
    real_t norm_squared = static_cast<real_t>(0);
    for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
      norm_squared += Q(row, col) * Q(row, col);
    }
    real_t norm = std::sqrt(static_cast<real_t>(norm_squared));
    if (checkRank && isFuzzyEqual(norm, 0.0)) {
      throw MathLibError("Matrix is rank deficient - cannot compute QR decomposition");
    }
    // R[col,col] = norm
    R(col, col) = static_cast<real_t>(norm);
    // Normalize: q_col = a_col / norm
    // Use multiplication by inverse for better performance
    if (!isFuzzyEqual(norm, 0.0)) {
      real_t norm_inv = static_cast<real_t>(1.0) / static_cast<real_t>(norm);
      for (index_t row = 1; row <= static_cast<index_t>(m); ++row) {
        Q(row, col) *= norm_inv;
      }
    }
  }
  QRDecomposition result;
  result.Q = Q;
  result.R = R;
  return result;
}


struct EigenDecomposition {
  Matrix<real_t> V; // eigenvectors
  Matrix<real_t> D; // eigenvalues on diagonal
};

struct EigenResult {
  Vector<real_t> values; // eigenvalues
  std::vector<Vector<real_t>> vectors; // eigenvectors
};

template<typename data_t>
Vector<real_t> eigvals(const Matrix<data_t>& mat, size_t max_iter = 1000, real_t tol = EPSILON)
{
  if (mat.rows() != mat.cols()) {
    throw MathLibError("Eigenvalue computation requires square matrix");
  }
  
  size_t n = mat.rows();
  if (n == 0) {
    return Vector<real_t>(0);
  }
  
  // QR algorithm: iteratively compute A_k = Q_k * R_k, then A_{k+1} = R_k * Q_k
  // The diagonal converges to eigenvalues
  Matrix<real_t> A(mat);
  
  for (size_t iter = 0; iter < max_iter; ++iter) {
    // Perform QR decomposition
    QRDecomposition qr_result = qr(A, false);
    
    // Form A_{k+1} = R * Q
    A = qr_result.R * qr_result.Q;
    
    // Check convergence: off-diagonal elements should approach zero
    real_t off_diag_sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
      for (size_t col = 0; col < n; ++col) {
        if (i != col) {
          off_diag_sum += std::abs(A[col * n + i]);
        }
      }
    }
    
    if (off_diag_sum < tol) {
      break;
    }
  }
  
  // Extract eigenvalues from diagonal
  Vector<real_t> eigenvalues(n);
  for (size_t i = 0; i < n; ++i) {
    eigenvalues[i] = A[i * n + i];
  }
  
  return eigenvalues;
}

template<typename data_t>
EigenDecomposition eigdecomp(const Matrix<data_t>& mat, size_t max_iter = 1000, real_t tol = EPSILON)
{
  if (mat.rows() != mat.cols()) {
    throw MathLibError("Eigenvalue computation requires square matrix");
  }
  
  size_t n = mat.rows();
  if (n == 0) {
    EigenDecomposition result;
    result.V = Matrix<real_t>(0, 0);
    result.D = Matrix<real_t>(0, 0);
    return result;
  }
  
  // QR algorithm with eigenvector accumulation
  Matrix<real_t> A(mat);
  Matrix<real_t> V = eye<real_t>(n); // Accumulate transformations
  
  for (size_t iter = 0; iter < max_iter; ++iter) {
    // Perform QR decomposition
    QRDecomposition qr_result = qr(A, false);
    
    // Accumulate eigenvectors: V = V * Q
    V = V * qr_result.Q;
    
    // Form A_{k+1} = R * Q
    A = qr_result.R * qr_result.Q;
    
    // Check convergence
    real_t off_diag_sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
      for (size_t col = 0; col < n; ++col) {
        if (i != col) {
          off_diag_sum += std::abs(A[col * n + i]);
        }
      }
    }
    
    if (off_diag_sum < tol) {
      break;
    }
  }
  
  EigenDecomposition result;
  result.V = V;
  result.D = A; // Diagonal matrix with eigenvalues
  
  return result;
}

template<typename data_t>
EigenResult eig(const Matrix<data_t>& mat, size_t max_iter = 1000, real_t tol = EPSILON)
{
  if (mat.rows() != mat.cols()) {
    throw MathLibError("Eigenvalue computation requires square matrix");
  }
  
  size_t n = mat.rows();
  EigenResult result;
  
  if (n == 0) {
    result.values = Vector<real_t>(0);
    result.vectors = std::vector<Vector<real_t>>();
    return result;
  }
  
  // Use eigdecomp to get both eigenvalues and eigenvectors
  EigenDecomposition decomp = eigdecomp(mat, max_iter, tol);
  
  // Extract eigenvalues from diagonal of D
  result.values = Vector<real_t>(n);
  for (size_t i = 0; i < n; ++i) {
    result.values[i] = decomp.D[i * n + i];
  }
  
  // Extract eigenvectors as columns of V
  result.vectors.resize(n);
  for (size_t col = 0; col < n; ++col) {
    result.vectors[col] = Vector<real_t>(n);
    for (size_t row = 0; row < n; ++row) {
      result.vectors[col][row] = decomp.V[col * n + row];
    }
  }
  
  return result;
}

inline Matrix<bool>
isFuzzyEqual(const Matrix<real_t>& lhs, const Matrix<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzyEqual(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Matrix<bool> isFuzzyGreater(const Matrix<real_t>& lhs, const Matrix<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzyGreater(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Matrix<bool> isFuzzySmaller(const Matrix<real_t>& lhs, const Matrix<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = isFuzzySmaller(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Matrix<bool> isStrictFuzzyGreater(const Matrix<real_t>& lhs, const Matrix<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = isStrictFuzzyGreater(lhs[i], rhs[i], epsilon);
  }
  return result;
}

inline Matrix<bool> isStrictFuzzySmaller(const Matrix<real_t>& lhs, const Matrix<real_t>& rhs, real_t epsilon = EPSILON)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw IncompatibleSizeError(lhs.rows(), lhs.cols(), rhs.rows(), rhs.cols());
  }
  Matrix<bool> result(lhs.rows(), lhs.cols());
  size_t len = lhs.length();
  for (size_t i = 0; i < len; ++i) {
    result[i] = isStrictFuzzySmaller(lhs[i], rhs[i], epsilon);
  }
  return result;
}

#pragma endregion

#ifdef __PRINT__
template <typename data_t>
std::string to_string(const Matrix<data_t>& mat, const size_t& decimal = 3, RoundingMethod_t method = NEAREST)
{
  if (mat.rows() == 0 || mat.cols() == 0) {
    return "[]";
  }

  // First pass: calculate maximum width needed for each column
  std::vector<size_t> col_widths(mat.cols(), 0);
  for (size_t col = 0; col < mat.cols(); ++col) {
    for (size_t row = 0; row < mat.rows(); ++row) {
      data_t rounded = mathlib::round(mat[col * mat.rows() + row], decimal, method);
      std::ostringstream temp;
      temp << std::fixed << std::setprecision(static_cast<int>(decimal)) << rounded;
      col_widths[col] = std::max(col_widths[col], temp.str().length());
    }
  }

  // Second pass: format with alignment
  std::ostringstream oss;
  oss << "[";
  for (size_t row = 0; row < mat.rows(); ++row) {
    if (row > 0)
      oss << " ";
    for (size_t col = 0; col < mat.cols(); ++col) {
      data_t rounded = mathlib::round(mat[col * mat.rows() + row], decimal, method);
      std::ostringstream temp;
      temp << std::fixed << std::setprecision(static_cast<int>(decimal)) << rounded;
      std::string value_str = temp.str();

      // Right-align the value within its column width
      oss << std::setw(static_cast<int>(col_widths[col])) << value_str;

      if (col < mat.cols() - 1)
        oss << ", ";
    }
    if (row < mat.rows() - 1)
      oss << ";\n";
  }
  oss << "]";
  return oss.str();
}

template <typename data_t>
std::ostream& operator<<(std::ostream& os, const Matrix<data_t>& mat)
{
  os << to_string(mat);
  return os;
}

template <typename data_t>
void print(const Matrix<data_t>& mat)
{
  std::cout << mat << '\n';
}
#endif

}  // namespace mathlib

#ifdef __PRINT__
namespace std
{
template <typename data_t>
std::string to_string(const mathlib::Matrix<data_t>& mat)
{
  return mathlib::to_string(mat);
}
}  // namespace std
#endif
