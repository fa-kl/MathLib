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
*/

namespace mathlib
{

#pragma region class Vector

template <typename data_t = real_t>
class Vector;

template <typename data_t>
class Matrix;

/// @brief A class for n-dimensional vectors
template <typename data_t>
class Vector
{
protected:
  /// @brief The underlying data array
  std::unique_ptr<data_t[]> m_data;
  /// @brief The size of the vector
  dim_t m_size;

public:
  /// @brief Creates an empty vector
  Vector() : m_data(nullptr), m_size(0, 0) {}

  /// @brief Creates a vector with n elements initialized to zero
  /// @param n Number of elements
  explicit Vector(size_t n) : m_data(std::make_unique<data_t[]>(n)), m_size(n, 1) {}

  /// @brief Creates a vector from an initializer list
  /// @param values Initializer list of values
  Vector(std::initializer_list<data_t> values)
      : m_data(std::make_unique<data_t[]>(values.size())), m_size(values.size(), 1)
  {
    std::copy(values.begin(), values.end(), m_data.get());
  }

  /// @brief Creates a vector from a matrix by flattening it
  /// @param mat Matrix to convert to vector
  Vector(const Matrix<data_t>& mat) : m_data(std::make_unique<data_t[]>(mat.length())), m_size(mat.length(), 1)
  {
    for (size_t i = 0; i < mat.length(); ++i) {
      m_data[i] = mat[i];
    }
  }

  /// @brief Copy constructor
  /// @param other Vector to copy
  Vector(const Vector& other)
      : m_data(other.length() > 0 ? std::make_unique<data_t[]>(other.length()) : nullptr), m_size(other.m_size)
  {
    if (m_data) {
      std::copy(other.m_data.get(), other.m_data.get() + length(), m_data.get());
    }
  }

  /// @brief Move constructor
  /// @param other Vector to move
  Vector(Vector&& other) noexcept : m_data(std::move(other.m_data)), m_size(other.m_size) { other.m_size = {0, 0}; }

  /// @brief Destructor
  ~Vector() = default;

  /// @brief Copy assignment operator
  /// @param other Vector to copy
  /// @returns Reference to this vector
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

  /// @brief Move assignment operator
  /// @param other Vector to move
  /// @returns Reference to this vector
  Vector& operator=(Vector&& other) noexcept
  {
    if (this != &other) {
      m_data = std::move(other.m_data);
      m_size = other.m_size;
      other.m_size = {0, 0};
    }
    return *this;
  }

  /// @brief Get the length (number of elements) of the vector
  /// @returns Length of the vector
  size_t length() const { return m_size.rows * m_size.cols; }

  /// @brief Get the size dimensions of the vector
  /// @returns Size as dim_t structure
  dim_t size() const { return m_size; }

  /// @brief Get the number of rows
  /// @returns Number of rows
  size_t rows() const { return m_size.rows; }

  /// @brief Get the number of columns
  /// @returns Number of columns
  size_t cols() const { return m_size.cols; }

  /// @brief Check if the vector is empty
  /// @returns True if empty, false otherwise
  bool isEmpty() const { return length() == 0; }

  /// @brief Access element at index i (0-based indexing)
  /// @param i Index
  /// @returns Reference to element at index i
  data_t& operator[](const size_t i)
  {
    if (i < length()) {
      return m_data[i];
    }
    throw IndexOutOfRangeError(i);
  }

  /// @brief Access element at index i (0-based indexing)
  /// @param i Index
  /// @returns Const reference to element at index i
  const data_t& operator[](const size_t i) const
  {
    if (i < length()) {
      return m_data[i];
    }
    throw IndexOutOfRangeError(i);
  }

  /// @brief Access element at index i (1-based indexing, supports negative indices)
  /// @param i Index (1-based, negative indices count from end)
  /// @returns Reference to element at index i
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

  /// @brief Access element at index i (1-based indexing, supports negative indices)
  /// @param i Index (1-based, negative indices count from end)
  /// @returns Const reference to element at index i
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

  /// @brief Unary plus operator
  /// @returns A copy of the vector
  Vector operator+() const { return Vector(*this); }

  /// @brief Unary minus operator (negate vector)
  /// @returns A negated copy of the vector
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

/// @brief Create a vector of zeros
/// @param n Length of the vector
/// @returns Vector filled with zeros
template <typename data_t>
Vector<data_t> zeros(size_t n)
{
  return Vector<data_t>(n);
}

/// @brief Create a vector of ones
/// @param n Length of the vector
/// @returns Vector filled with ones
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

/// @brief Add two vectors element-wise
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result of element-wise addition
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

/// @brief Add two vectors of different types element-wise
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result of element-wise addition with appropriate result type
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

/// @brief Add a scalar to all elements of a vector
/// @param lhs Vector
/// @param rhs Scalar value
/// @returns Result vector
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

/// @brief Add a scalar to all elements of a vector
/// @param lhs Scalar value
/// @param rhs Vector
/// @returns Result vector
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

/// @brief Subtract two vectors element-wise
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result of element-wise subtraction
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

/// @brief Subtract two vectors of different types element-wise
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result of element-wise subtraction with appropriate result type
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

/// @brief Subtract a scalar from all elements of a vector
/// @param lhs Vector
/// @param rhs Scalar value
/// @returns Result vector
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

/// @brief Subtract a vector from a scalar
/// @param lhs Scalar value
/// @param rhs Vector
/// @returns Result vector
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

/// @brief Compute dot product of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Dot product (scalar value)
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

/// @brief Compute dot product of two vectors of different types
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Dot product (scalar value)
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

/// @brief Multiply all elements of a vector by a scalar
/// @param lhs Vector
/// @param rhs Scalar value
/// @returns Result vector
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

/// @brief Multiply all elements of a vector by a scalar
/// @param lhs Scalar value
/// @param rhs Vector
/// @returns Result vector
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

/// @brief Divide all elements of a vector by a scalar
/// @param lhs Vector
/// @param rhs Scalar value (divisor)
/// @returns Result vector
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

/// @brief Divide a scalar by all elements of a vector
/// @param lhs Scalar value (dividend)
/// @param rhs Vector (divisor)
/// @returns Result vector
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

/// @brief Element-wise logical AND of two boolean vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result boolean vector
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

/// @brief Element-wise logical OR of two boolean vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result boolean vector
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

/// @brief Element-wise logical XOR of two boolean vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Result boolean vector
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

/// @brief Element-wise equality comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise inequality comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise greater-than comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise greater-than-or-equal comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise less-than comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise less-than-or-equal comparison of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Boolean vector with comparison results
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

/// @brief Get the length of a vector
/// @param vec Input vector
/// @returns Length of the vector
template <typename data_t>
size_t length(const Vector<data_t>& vec)
{
  return vec.length();
}

/// @brief Get the size dimensions of a vector
/// @param vec Input vector
/// @returns Size as dim_t structure
template <typename data_t>
dim_t size(const Vector<data_t>& vec)
{
  return vec.size();
}

/// @brief Check if a vector is empty
/// @param vec Input vector
/// @returns True if empty, false otherwise
template <typename data_t>
bool isEmpty(const Vector<data_t>& vec)
{
  return vec.isEmpty();
}

/// @brief Transpose a vector (convert column vector to row vector or vice versa)
/// @param vec Input vector
/// @returns Transposed matrix
template <typename data_t>
Matrix<data_t> transpose(const Vector<data_t>& vec)
{
  Matrix<data_t> result(vec.cols(), vec.rows());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = vec[i];
  }
  return result;
}

/// @brief Horizontally concatenate two vectors into a matrix
/// @param lhs Left vector (first column)
/// @param rhs Right vector (second column)
/// @returns Matrix with two columns
template <typename data_t>
Matrix<data_t> hcat(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError("Vectors must have the same length for horizontal concatenation");
  }
  size_t rows = lhs.length();
  Matrix<data_t> result(rows, 2);

  // Copy lhs into first column using 1-based indexing
  for (index_t i = 1; i <= static_cast<index_t>(rows); ++i) {
    result(i, 1) = lhs(i);
  }

  // Copy rhs into second column using 1-based indexing
  for (index_t i = 1; i <= static_cast<index_t>(rows); ++i) {
    result(i, 2) = rhs(i);
  }

  return result;
}

/// @brief Vertically concatenate two vectors
/// @param lhs First vector
/// @param rhs Second vector
/// @returns Concatenated vector
template <typename data_t>
Vector<data_t> vcat(const Vector<data_t>& lhs, const Vector<data_t>& rhs)
{
  size_t lhs_len = lhs.length();
  size_t rhs_len = rhs.length();
  Vector<data_t> result(lhs_len + rhs_len);

  // Copy lhs elements using 1-based indexing
  for (index_t i = 1; i <= static_cast<index_t>(lhs_len); ++i) {
    result(i) = lhs(i);
  }

  // Copy rhs elements using 1-based indexing
  for (index_t i = 1; i <= static_cast<index_t>(rhs_len); ++i) {
    result(static_cast<index_t>(lhs_len) + i) = rhs(i);
  }

  return result;
}

/// @brief Compute the p-norm of a vector
/// @param vec Input vector
/// @param p Norm parameter (default: 2 for Euclidean norm)
/// @returns p-norm of the vector
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

/// @brief Normalize a vector (divide by its norm)
/// @param vec Input vector
/// @returns Normalized vector
template <typename data_t>
Vector<data_t> normalize(const Vector<data_t>& vec)
{
  real_t n = norm(vec);
  if (isFuzzyEqual(n, 0.0)) {
    throw MathLibError("Cannot normalize a zero-length vector.");
  }
  return vec / n;
}

/// @brief Compute absolute value of each element in a vector
/// @param vec Input vector
/// @returns Vector with absolute values
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

/// @brief Compute exponential of each element in a vector
/// @param vec Input vector
/// @returns Vector with exponential values
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

/// @brief Compute natural logarithm of each element in a vector
/// @param vec Input vector
/// @returns Vector with logarithmic values
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

/// @brief Compute sine of each element in a vector
/// @param vec Input vector
/// @returns Vector with sine values
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

/// @brief Compute cosine of each element in a vector
/// @param vec Input vector
/// @returns Vector with cosine values
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

/// @brief Compute tangent of each element in a vector
/// @param vec Input vector
/// @returns Vector with tangent values
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

/// @brief Compute arcsine of each element in a vector
/// @param vec Input vector
/// @returns Vector with arcsine values
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

/// @brief Compute arccosine of each element in a vector
/// @param vec Input vector
/// @returns Vector with arccosine values
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

/// @brief Compute arctangent of each element in a vector
/// @param vec Input vector
/// @returns Vector with arctangent values
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

/// @brief Raise each element of a vector to a power
/// @param vec Input vector
/// @param exponent Power exponent
/// @returns Vector with powered values
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

/// @brief Compute square root of each element in a vector
/// @param vec Input vector
/// @returns Vector with square root values
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

/// @brief Round each element of a vector to a specified decimal place
/// @param vec Input vector
/// @param decimal Number of decimal places
/// @param method Rounding method (default: NEAREST)
/// @returns Vector with rounded values
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

/// @brief Compute sum of all elements in a vector
/// @param vec Input vector
/// @returns Sum of elements
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

/// @brief Element-wise sum of two vectors of different types
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Vector with element-wise sum
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

/// @brief Compute differences between consecutive elements
/// @param vec Input vector
/// @returns Vector of differences (length n-1)
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

/// @brief Element-wise difference of two vectors of different types
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Vector with element-wise difference
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

/// @brief Element-wise multiplication of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Vector with element-wise products
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

/// @brief Element-wise division of two vectors
/// @param lhs Left-hand side vector (dividend)
/// @param rhs Right-hand side vector (divisor)
/// @returns Vector with element-wise quotients
template <typename T1, typename T2>
auto eldiv(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (lhs.length() != rhs.length()) {
    throw IncompatibleSizeError();
  }
  using result_t = std::common_type_t<T1, T2>;
  size_t len = lhs.length();
  Vector<result_t> result(len);
  for (size_t i = 0; i < len; ++i) {
    result[i] = static_cast<result_t>(lhs[i]) / static_cast<result_t>(rhs[i]);
  }
  return result;
}

/// @brief Compute dot product of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Dot product (scalar value)
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

/// @brief Compute cross product of two 3D vectors
/// @param lhs Left-hand side vector (must be 3D)
/// @param rhs Right-hand side vector (must be 3D)
/// @returns Cross product vector
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

/// @brief Compute outer product of two vectors
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @returns Matrix representing the outer product
template <typename T1, typename T2>
auto outer(const Vector<T1>& lhs, const Vector<T2>& rhs) -> Matrix<decltype(std::declval<T1>() + std::declval<T2>())>
{
  using result_t = std::common_type_t<T1, T2>;
  const size_t m = lhs.length();
  const size_t n = rhs.length();

  Matrix<result_t> result(m, n);

  for (size_t i = 0; i < m; ++i) {
    for (size_t col = 0; col < n; ++col) {
      result(static_cast<index_t>(i + 1), static_cast<index_t>(col + 1)) =
          static_cast<result_t>(lhs[i]) * static_cast<result_t>(rhs[col]);
    }
  }

  return result;
}

/// @brief Find maximum element in a vector
/// @param vec Input vector
/// @returns Maximum value
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

/// @brief Find minimum element in a vector
/// @param vec Input vector
/// @returns Minimum value
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

/// @brief Compute mean (average) of vector elements
/// @param vec Input vector
/// @returns Mean value
template <typename data_t>
real_t mean(const Vector<data_t>& vec)
{
  return sum(vec) / static_cast<real_t>(vec.length());
}

/// @brief Compute standard deviation of vector elements
/// @param vec Input vector
/// @returns Standard deviation
template <typename data_t>
real_t std(const Vector<data_t>& vec)
{
  return std::sqrt(var(vec));
}

/// @brief Compute variance of vector elements
/// @param vec Input vector
/// @returns Variance
template <typename data_t>
real_t var(const Vector<data_t>& vec)
{
  real_t factor = static_cast<real_t>(1.0) / (static_cast<real_t>(vec.length()) - static_cast<real_t>(1.0));
  Vector<real_t> v = abs(vec - mean(vec));
  return factor * sum(elmul(v, v));
}

/// @brief Check if all elements in a boolean vector are true
/// @param vec Boolean vector
/// @returns True if all elements are true, false otherwise
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

/// @brief Check if any element in a boolean vector is true
/// @param vec Boolean vector
/// @returns True if at least one element is true, false otherwise
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

/// @brief Element-wise fuzzy equality comparison with tolerance
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @param epsilon Tolerance value (default: EPSILON)
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise strict fuzzy greater-than comparison with tolerance
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @param epsilon Tolerance value (default: EPSILON)
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise fuzzy greater-than-or-equal comparison with tolerance
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @param epsilon Tolerance value (default: EPSILON)
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise strict fuzzy less-than comparison with tolerance
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @param epsilon Tolerance value (default: EPSILON)
/// @returns Boolean vector with comparison results
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

/// @brief Element-wise fuzzy less-than-or-equal comparison with tolerance
/// @param lhs Left-hand side vector
/// @param rhs Right-hand side vector
/// @param epsilon Tolerance value (default: EPSILON)
/// @returns Boolean vector with comparison results
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
/// @brief Convert a vector to a string representation
/// @param vec Input vector
/// @param decimal Number of decimal places (default: 3)
/// @param method Rounding method (default: NEAREST)
/// @returns String representation of the vector
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

/// @brief Output stream operator for vectors
/// @param os Output stream
/// @param vec Vector to output
/// @returns Output stream
template <typename data_t>
std::ostream& operator<<(std::ostream& os, const Vector<data_t>& vec)
{
  os << to_string(vec);
  return os;
}

/// @brief Print a vector to standard output
/// @param vec Vector to print
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
/// @brief Convert a vector to a string (std namespace overload)
/// @param vec Input vector
/// @returns String representation of the vector
template <typename data_t>
std::string to_string(const mathlib::Vector<data_t>& vec)
{
  return mathlib::to_string(vec);
}

/// @brief Create a std::vector from this Vector
template <typename data_t>
std::vector<data_t> to_vector(const mathlib::Vector<data_t>& vec)
{
  size_t len = vec.length();
  ;
  std::vector<data_t> result;
  result.reserve(len);
  for (size_t i = 0; i < len; ++i) {
    result.push_back(vec[i]);
  }
  return result;
}

#endif

}  // namespace std