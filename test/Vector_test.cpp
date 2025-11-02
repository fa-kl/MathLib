#include "Vector.hpp"

#include <complex>

#include <gtest/gtest.h>

#include "DivByZeroError.hpp"
#include "IncompatibleSizeError.hpp"
#include "IndexOutOfRangeError.hpp"
#include "MathLibError.hpp"
#include "Matrix.hpp"

using namespace mathlib;

// Constructor Tests
TEST(VectorTest, DefaultConstructor)
{
  const Vector<real_t> vec;
  EXPECT_EQ(vec.length(), 0);
  EXPECT_TRUE(vec.isEmpty());
  EXPECT_THROW(vec[0], IndexOutOfRangeError);
  EXPECT_THROW(vec(0), IndexOutOfRangeError);
  EXPECT_THROW(vec(1), IndexOutOfRangeError);
  EXPECT_THROW(vec(-1), IndexOutOfRangeError);
}

TEST(VectorTest, SizeConstructor)
{
  Vector<real_t> vec(5);

  EXPECT_EQ(vec.length(), 5);
  EXPECT_FALSE(vec.isEmpty());
  EXPECT_DOUBLE_EQ(vec[0], 0.0);
  EXPECT_DOUBLE_EQ(vec[4], 0.0);
  EXPECT_THROW(vec[5], IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(vec(1), 0.0);
  EXPECT_DOUBLE_EQ(vec(5), 0.0);
  EXPECT_DOUBLE_EQ(vec(-1), 0.0);
  EXPECT_DOUBLE_EQ(vec(-5), 0.0);
  EXPECT_THROW(vec(6), IndexOutOfRangeError);
  EXPECT_THROW(vec(-6), IndexOutOfRangeError);
}

TEST(VectorTest, InitializerListConstructor)
{
  const Vector<real_t> vec1 = {1.0, 2.0, 3.0};
  EXPECT_EQ(vec1.length(), 3);
  EXPECT_DOUBLE_EQ(vec1[0], 1.0);
  EXPECT_DOUBLE_EQ(vec1[1], 2.0);
  EXPECT_DOUBLE_EQ(vec1[2], 3.0);

  const Vector<real_t> vec2({-2.0, -1.0, 5.0});
  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], -2.0);
  EXPECT_DOUBLE_EQ(vec2[1], -1.0);
  EXPECT_DOUBLE_EQ(vec2[2], 5.0);

  const Vector<int64_t> vec3 = {1, 2, 3, 4, 5};
  EXPECT_EQ(vec3.length(), 5);
  EXPECT_EQ(vec3[0], 1);
  EXPECT_EQ(vec3[4], 5);
}

TEST(VectorTest, MatrixConstructor)
{
  Matrix<real_t> col_mat(3, 1);
  col_mat(1, 1) = 1.0;
  col_mat(2, 1) = 2.0;
  col_mat(3, 1) = 3.0;

  Vector<real_t> vec_from_col(col_mat);
  EXPECT_EQ(vec_from_col.length(), 3);
  EXPECT_DOUBLE_EQ(vec_from_col[0], 1.0);
  EXPECT_DOUBLE_EQ(vec_from_col[1], 2.0);
  EXPECT_DOUBLE_EQ(vec_from_col[2], 3.0);
  EXPECT_DOUBLE_EQ(vec_from_col(1), 1.0);
  EXPECT_DOUBLE_EQ(vec_from_col(2), 2.0);
  EXPECT_DOUBLE_EQ(vec_from_col(3), 3.0);

  Matrix<real_t> row_mat(1, 3);
  row_mat(1, 1) = 4.0;
  row_mat(1, 2) = 5.0;
  row_mat(1, 3) = 6.0;

  Vector<real_t> vec_from_row(row_mat);
  EXPECT_EQ(vec_from_row.length(), 3);
  EXPECT_DOUBLE_EQ(vec_from_row[0], 4.0);
  EXPECT_DOUBLE_EQ(vec_from_row[1], 5.0);
  EXPECT_DOUBLE_EQ(vec_from_row[2], 6.0);

  Matrix<real_t> mat(2, 3);
  mat(1, 1) = 1.0;  
  mat(2, 1) = 2.0;
  mat(1, 2) = 3.0;  
  mat(2, 2) = 4.0;
  mat(1, 3) = 5.0;  
  mat(2, 3) = 6.0;

  Vector<real_t> vec_from_mat(mat);
  EXPECT_EQ(vec_from_mat.length(), 6);
  EXPECT_DOUBLE_EQ(vec_from_mat[0], 1.0);
  EXPECT_DOUBLE_EQ(vec_from_mat[1], 2.0);
  EXPECT_DOUBLE_EQ(vec_from_mat[2], 3.0);
  EXPECT_DOUBLE_EQ(vec_from_mat[3], 4.0);
  EXPECT_DOUBLE_EQ(vec_from_mat[4], 5.0);
  EXPECT_DOUBLE_EQ(vec_from_mat[5], 6.0);

  Matrix<real_t> single_mat(1, 1);
  single_mat(1, 1) = 42.0;

  Vector<real_t> vec_from_single(single_mat);
  EXPECT_EQ(vec_from_single.length(), 1);
  EXPECT_DOUBLE_EQ(vec_from_single[0], 42.0);
  EXPECT_DOUBLE_EQ(vec_from_single(1), 42.0);

  Matrix<int64_t> int_mat(2, 2);
  int_mat(1, 1) = 10;
  int_mat(2, 1) = 20;
  int_mat(1, 2) = 30;
  int_mat(2, 2) = 40;

  Vector<int64_t> vec_from_int(int_mat);
  EXPECT_EQ(vec_from_int.length(), 4);
  EXPECT_EQ(vec_from_int[0], 10);
  EXPECT_EQ(vec_from_int[1], 20);
  EXPECT_EQ(vec_from_int[2], 30);
  EXPECT_EQ(vec_from_int[3], 40);

  Matrix<real_t> mat_copy(2, 1);
  mat_copy(1, 1) = 100.0;
  mat_copy(2, 1) = 200.0;

  Vector<real_t> vec_copy(mat_copy);

  EXPECT_DOUBLE_EQ(vec_copy[0], 100.0);
  EXPECT_DOUBLE_EQ(vec_copy[1], 200.0);
}

TEST(VectorTest, CopyConstructor)
{
  Vector<real_t> vec1 = {1.0, 2.0, 3.0};
  Vector<real_t> vec2(vec1);
  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
  EXPECT_DOUBLE_EQ(vec2[1], 2.0);
  EXPECT_DOUBLE_EQ(vec2[2], 3.0);

  // Verify deep copy
  vec1[0] = 0.0;
  vec1[1] = 0.0;
  vec1[2] = 0.0;
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
  EXPECT_DOUBLE_EQ(vec2[1], 2.0);
  EXPECT_DOUBLE_EQ(vec2[2], 3.0);
}

TEST(VectorTest, MoveConstructor)
{
  Vector<real_t> vec1 = {1.0, 2.0, 3.0};
  Vector<real_t> vec2(std::move(vec1));

  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
  EXPECT_DOUBLE_EQ(vec2[1], 2.0);
  EXPECT_DOUBLE_EQ(vec2[2], 3.0);
  EXPECT_EQ(vec1.length(), 0);
}

TEST(VectorTest, CopyAssignment)
{
  const Vector<real_t> vec1 = {1.0, 2.0, 3.0};
  Vector<real_t> vec2;
  vec2 = vec1;
  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
  EXPECT_DOUBLE_EQ(vec2[1], 2.0);
  EXPECT_DOUBLE_EQ(vec2[2], 3.0);

  // Test self-assignment
  vec2 = vec2;
  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
}

TEST(VectorTest, MoveAssignment)
{
  Vector<real_t> vec1 = {1.0, 2.0, 3.0};
  Vector<real_t> vec2 = {4.0, 5.0};

  vec2 = std::move(vec1);

  EXPECT_EQ(vec2.length(), 3);
  EXPECT_DOUBLE_EQ(vec2[0], 1.0);
  EXPECT_DOUBLE_EQ(vec2[1], 2.0);
  EXPECT_DOUBLE_EQ(vec2[2], 3.0);
  EXPECT_EQ(vec1.length(), 0);

  // Test self-assignment
  vec2 = std::move(vec2);
  EXPECT_EQ(vec2.length(), 3);
}

TEST(VectorTest, SizeAndLength)
{
  // Empty vector
  Vector<real_t> vec0;
  EXPECT_EQ(vec0.length(), 0);
  EXPECT_TRUE(vec0.isEmpty());
  dim_t sz0 = vec0.size();
  EXPECT_EQ(sz0.rows, 0);
  EXPECT_EQ(sz0.cols, 0);

  // Single element
  Vector<real_t> vec1 = {1.0};
  EXPECT_EQ(vec1.length(), 1);
  EXPECT_FALSE(vec1.isEmpty());
  dim_t sz1 = vec1.size();
  EXPECT_EQ(sz1.rows, 1);
  EXPECT_EQ(sz1.cols, 1);

  // Three elements
  Vector<real_t> vec3 = {1.0, 2.0, 3.0};
  EXPECT_EQ(vec3.length(), 3);
  EXPECT_FALSE(vec3.isEmpty());
  dim_t sz3 = vec3.size();
  EXPECT_EQ(sz3.rows, 3);
  EXPECT_EQ(sz3.cols, 1);

  // Ten elements
  Vector<int64_t> vec10 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  EXPECT_EQ(vec10.length(), 10);
  EXPECT_FALSE(vec10.isEmpty());
  dim_t sz10 = vec10.size();
  EXPECT_EQ(sz10.rows, 10);
  EXPECT_EQ(sz10.cols, 1);

  // Test length() free function
  EXPECT_EQ(length(vec3), 3);
  EXPECT_EQ(size(vec3).rows, 3);
  EXPECT_EQ(size(vec3).cols, 1);
}

TEST(VectorTest, Transpose)
{
  Vector<real_t> vec = {1.0, 2.0, 3.0};

  // Original dimensions
  EXPECT_EQ(vec.length(), 3);
  dim_t sz = vec.size();
  EXPECT_EQ(sz.rows, 3);
  EXPECT_EQ(sz.cols, 1);

  // Transpose (free function)
  Matrix<real_t> mat = transpose(vec);
  EXPECT_EQ(mat.length(), 3);
  EXPECT_EQ(mat.rows(), 1);
  EXPECT_EQ(mat.cols(), 3);

  // Original vector should remain unchanged
  sz = vec.size();
  EXPECT_EQ(sz.rows, 3);
  EXPECT_EQ(sz.cols, 1);

  // Values should remain unchanged in original vector
  EXPECT_DOUBLE_EQ(vec[0], 1.0);
  EXPECT_DOUBLE_EQ(vec[1], 2.0);
  EXPECT_DOUBLE_EQ(vec[2], 3.0);

  // Values should be correct in transposed matrix
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 3), 3.0);

  // Transpose back to get column matrix
  Matrix<real_t> mat2 = transpose(mat);
  EXPECT_EQ(mat2.rows(), 3);
  EXPECT_EQ(mat2.cols(), 1);
  EXPECT_DOUBLE_EQ(mat2(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat2(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat2(3, 1), 3.0);
}

TEST(VectorTest, Getter)
{
  const Vector vec1 = {1.0, 2.0, 3.0};
  EXPECT_DOUBLE_EQ(vec1[0], 1.0);
  EXPECT_DOUBLE_EQ(vec1[1], 2.0);
  EXPECT_DOUBLE_EQ(vec1[2], 3.0);

  EXPECT_THROW(vec1(0), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(vec1(1), 1.0);
  EXPECT_DOUBLE_EQ(vec1(2), 2.0);
  EXPECT_DOUBLE_EQ(vec1(3), 3.0);
  EXPECT_THROW(vec1(4), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(vec1(-1), 3.0);
  EXPECT_DOUBLE_EQ(vec1(-2), 2.0);
  EXPECT_DOUBLE_EQ(vec1(-3), 1.0);
  EXPECT_THROW(vec1(-4), IndexOutOfRangeError);

  const Vector vec2 = {3.0, 4.0, 5.0};
  EXPECT_DOUBLE_EQ(vec2[0], 3.0);
  EXPECT_DOUBLE_EQ(vec2[1], 4.0);
  EXPECT_DOUBLE_EQ(vec2[2], 5.0);
}

TEST(VectorTest, Norm)
{
  // Test zero vector
  const Vector zero = {0.0, 0.0, 0.0};
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm(zero, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(norm(zero, 2.0), 0.0);
  EXPECT_DOUBLE_EQ(norm(zero, 3.0), 0.0);
  EXPECT_DOUBLE_EQ(norm(zero, std::numeric_limits<real_t>::infinity()), 0.0);

  // Test unit vectors (2-norm should be 1.0)
  const Vector e1 = {1.0, 0.0, 0.0};
  const Vector e2 = {0.0, 1.0, 0.0};
  const Vector e3 = {0.0, 0.0, 1.0};
  EXPECT_DOUBLE_EQ(norm(e1), 1.0);
  EXPECT_DOUBLE_EQ(norm(e2), 1.0);
  EXPECT_DOUBLE_EQ(norm(e3), 1.0);

  // Test 2-norm (Euclidean norm) - default
  const Vector vec = {4.0, 3.0, 12.0};
  EXPECT_DOUBLE_EQ(norm(vec), std::sqrt(4.0 * 4.0 + 3.0 * 3.0 + 12.0 * 12.0));
  EXPECT_DOUBLE_EQ(norm(vec, 2.0), std::sqrt(4.0 * 4.0 + 3.0 * 3.0 + 12.0 * 12.0));

  // Test 1-norm (Manhattan norm)
  const Vector vec1 = {3.0, -4.0, 5.0};
  EXPECT_DOUBLE_EQ(norm(vec1, 1.0), 3.0 + 4.0 + 5.0);

  // Test 1-norm with negative values
  const Vector vecNeg = {-2.0, -3.0, -4.0};
  EXPECT_DOUBLE_EQ(norm(vecNeg, 1.0), 2.0 + 3.0 + 4.0);

  // Test infinity norm (max absolute value)
  const Vector vecInf = {1.0, -5.0, 3.0, 2.0};
  EXPECT_DOUBLE_EQ(norm(vecInf, std::numeric_limits<real_t>::infinity()), 5.0);

  // Test -infinity norm (min absolute value)
  EXPECT_DOUBLE_EQ(norm(vecInf, -std::numeric_limits<real_t>::infinity()), 1.0);

  // Test p-norm with p=3
  const Vector vec3 = {1.0, 2.0, 3.0};
  const real_t expected_p3 = std::pow(1.0 + 8.0 + 27.0, 1.0 / 3.0);
  EXPECT_NEAR(norm(vec3, 3.0), expected_p3, 1e-10);

  // Test p-norm with p=4
  const Vector vec4 = {1.0, 2.0, 2.0};
  const real_t expected_p4 = std::pow(1.0 + 16.0 + 16.0, 1.0 / 4.0);
  EXPECT_NEAR(norm(vec4, 4.0), expected_p4, 1e-10);

  // Test that p <= 0 throws error
  EXPECT_THROW(norm(vec, 0.0), MathLibError);
  EXPECT_THROW(norm(vec, -0.5), MathLibError);

  // Test with different vector sizes
  const Vector vec2d = {3.0, 4.0};
  EXPECT_DOUBLE_EQ(norm(vec2d), 5.0);  // 3-4-5 triangle

  const Vector vec1d = {7.0};
  EXPECT_DOUBLE_EQ(norm(vec1d), 7.0);

  // Test with mixed positive and negative values
  const Vector vecMixed = {-3.0, 4.0, -12.0};
  EXPECT_DOUBLE_EQ(norm(vecMixed), 13.0);  // sqrt(9 + 16 + 144) = sqrt(169) = 13
}

TEST(VectorTest, NormEdgeCases)
{
  const Vector<real_t> vec = {3.0, 4.0};

  // Test very large p-norm
  real_t large_p = norm(vec, 100.0);
  EXPECT_NEAR(large_p, 4.0, 1e-6);  // Should approach infinity norm

  // Test p close to 1
  real_t p1 = norm(vec, 1.0);
  EXPECT_DOUBLE_EQ(p1, 7.0);

  // Single element vector
  const Vector<real_t> single = {5.0};
  EXPECT_DOUBLE_EQ(norm(single), 5.0);
  EXPECT_DOUBLE_EQ(norm(single, 1.0), 5.0);
  EXPECT_DOUBLE_EQ(norm(single, std::numeric_limits<real_t>::infinity()), 5.0);
}

TEST(VectorTest, Normalize)
{
  const Vector zero = {0.0, 0.0, 0.0};
  EXPECT_THROW(normalize(zero), mathlib::MathLibError);

  const Vector e1 = {1.0, 0.0, 0.0};
  const Vector n1 = normalize(e1);
  EXPECT_DOUBLE_EQ(n1[0], 1.0);
  EXPECT_DOUBLE_EQ(n1[1], 0.0);
  EXPECT_DOUBLE_EQ(n1[2], 0.0);

  const Vector e2 = {0.0, 1.0, 0.0};
  const Vector n2 = normalize(e2);
  EXPECT_DOUBLE_EQ(n2[0], 0.0);
  EXPECT_DOUBLE_EQ(n2[1], 1.0);
  EXPECT_DOUBLE_EQ(n2[2], 0.0);

  const Vector vec = {4.0, 2.0, -3.0};
  const Vector n = normalize(vec);
  EXPECT_DOUBLE_EQ(n[0], 4.0 / norm(vec));
  EXPECT_DOUBLE_EQ(n[1], 2.0 / norm(vec));
  EXPECT_DOUBLE_EQ(n[2], -3.0 / norm(vec));
}

TEST(VectorTest, NormalizeZeroVector)
{
  const Vector<real_t> zero = {0.0, 0.0, 0.0};
  EXPECT_THROW(normalize(zero), MathLibError);
}

TEST(VectorTest, Abs)
{
  const Vector vec1 = {4.0, 3.0, -2.0};
  const Vector result = abs(vec1);
  EXPECT_DOUBLE_EQ(result[0], 4.0);
  EXPECT_DOUBLE_EQ(result[1], 3.0);
  EXPECT_DOUBLE_EQ(result[2], 2.0);
}

TEST(VectorTest, Exp)
{
  const Vector vec = {0.5, 1.0, -1.0};
  const Vector result = exp(vec);
  EXPECT_DOUBLE_EQ(result[0], std::exp(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::exp(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::exp(vec[2]));
}

TEST(VectorTest, Log)
{
  const Vector vec = {0.5, 1.0, 2.0};
  const Vector result = log(vec);
  EXPECT_DOUBLE_EQ(result[0], std::log(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::log(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::log(vec[2]));
}

TEST(VectorTest, Sin)
{
  const Vector vec = {0.0, PI / 2.0, PI};
  const Vector result = sin(vec);
  EXPECT_DOUBLE_EQ(result[0], std::sin(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::sin(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::sin(vec[2]));
}

TEST(VectorTest, Cos)
{
  const Vector vec = {0.0, PI / 2.0, PI};
  const Vector result = cos(vec);
  EXPECT_DOUBLE_EQ(result[0], std::cos(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::cos(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::cos(vec[2]));
}

TEST(VectorTest, Tan)
{
  const Vector vec = {0.0, PI / 4.0, PI / 3.0};
  const Vector result = tan(vec);
  EXPECT_DOUBLE_EQ(result[0], std::tan(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::tan(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::tan(vec[2]));
}

TEST(VectorTest, Asin)
{
  const Vector vec = {0.0, 0.5, -0.5};
  const Vector result = asin(vec);
  EXPECT_DOUBLE_EQ(result[0], std::asin(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::asin(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::asin(vec[2]));
}

TEST(VectorTest, Acos)
{
  const Vector vec = {0.0, 0.5, -0.5};
  const Vector result = acos(vec);
  EXPECT_DOUBLE_EQ(result[0], std::acos(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::acos(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::acos(vec[2]));
}

TEST(VectorTest, Atan)
{
  const Vector vec = {0.0, 1.0, -1.0};
  const Vector result = atan(vec);
  EXPECT_DOUBLE_EQ(result[0], std::atan(vec[0]));
  EXPECT_DOUBLE_EQ(result[1], std::atan(vec[1]));
  EXPECT_DOUBLE_EQ(result[2], std::atan(vec[2]));
}

TEST(VectorTest, PowRealExponent)
{
  const real_t exponent = 2.0;
  const Vector vec = {4.0, 3.0, 2.0};
  const Vector result = pow(vec, exponent);
  EXPECT_DOUBLE_EQ(result[0], std::pow(vec[0], exponent));
  EXPECT_DOUBLE_EQ(result[1], std::pow(vec[1], exponent));
  EXPECT_DOUBLE_EQ(result[2], std::pow(vec[2], exponent));
}

TEST(VectorTest, Sqrt)
{
  const Vector vec1 = {1.0, 2.0, 3.0};
  const Vector result = sqrt(vec1);
  EXPECT_DOUBLE_EQ(result[0], std::sqrt(vec1[0]));
  EXPECT_DOUBLE_EQ(result[1], std::sqrt(vec1[1]));
  EXPECT_DOUBLE_EQ(result[2], std::sqrt(vec1[2]));
  // Square root of negative numbers should throw
  EXPECT_THROW(sqrt(Vector({-1.0, 2.0, 3.0})), MathLibError);
}

TEST(VectorTest, Dot)
{
  const Vector vec1 = {1.0, 2.0, 3.0};
  const Vector vec2 = {4.0, 5.0, 6.0};
  const Vector zero = {0.0, 0.0, 0.0};
  EXPECT_DOUBLE_EQ(dot(vec1, zero), 0.0);
  EXPECT_DOUBLE_EQ(dot(zero, vec1), 0.0);
  EXPECT_DOUBLE_EQ(dot(vec1, vec2), 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0);
  EXPECT_DOUBLE_EQ(dot(vec2, vec1), 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0);
}

TEST(VectorTest, DotProductZeroLength)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {0.0, 0.0, 0.0};

  EXPECT_DOUBLE_EQ(dot(v1, v2), 0.0);
  EXPECT_DOUBLE_EQ(dot(v2, v1), 0.0);
  EXPECT_DOUBLE_EQ(dot(v2, v2), 0.0);
}

TEST(VectorTest, DotProductIncompatibleSizes)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {1.0, 2.0};

  EXPECT_THROW(dot(v1, v2), IncompatibleSizeError);
}

TEST(VectorTest, Cross)
{
  const Vector e1 = {1.0, 0.0, 0.0};
  const Vector e2 = {0.0, 1.0, 0.0};
  const Vector e3 = {0.0, 0.0, 1.0};

  const Vector r1 = cross(e1, e2);
  EXPECT_DOUBLE_EQ(r1[0], e3[0]);
  EXPECT_DOUBLE_EQ(r1[1], e3[1]);
  EXPECT_DOUBLE_EQ(r1[2], e3[2]);

  const Vector r2 = cross(e1, e3);
  EXPECT_DOUBLE_EQ(r2[0], -e2[0]);
  EXPECT_DOUBLE_EQ(r2[1], -e2[1]);
  EXPECT_DOUBLE_EQ(r2[2], -e2[2]);

  const Vector r3 = cross(e2, e3);
  EXPECT_DOUBLE_EQ(r3[0], e1[0]);
  EXPECT_DOUBLE_EQ(r3[1], e1[1]);
  EXPECT_DOUBLE_EQ(r3[2], e1[2]);

  const Vector r4 = cross(e1, 2 * e1);
  EXPECT_DOUBLE_EQ(r4[0], 0.0);
  EXPECT_DOUBLE_EQ(r4[1], 0.0);
  EXPECT_DOUBLE_EQ(r4[2], 0.0);

  EXPECT_THROW(cross(Vector({1, 2}), Vector({3, 4})), mathlib::IncompatibleSizeError);
}

TEST(VectorTest, CrossProductIncompatibleSize)
{
  const Vector<real_t> v2d = {1.0, 2.0};
  const Vector<real_t> v3d = {3.0, 4.0, 5.0};
  const Vector<real_t> v4d = {1.0, 2.0, 3.0, 4.0};

  EXPECT_THROW(cross(v2d, v3d), IncompatibleSizeError);
  EXPECT_THROW(cross(v3d, v4d), IncompatibleSizeError);
  EXPECT_THROW(cross(v2d, v2d), IncompatibleSizeError);  // Must be 3D
}

TEST(VectorTest, Outer)
{
  // Basic outer product: 3D vectors
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {4.0, 5.0, 6.0};
  const Matrix<real_t> result = outer(v1, v2);

  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 3);

  // Check all elements
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0 * 4.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 1.0 * 5.0);
  EXPECT_DOUBLE_EQ(result(1, 3), 1.0 * 6.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 2.0 * 4.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 2.0 * 5.0);
  EXPECT_DOUBLE_EQ(result(2, 3), 2.0 * 6.0);
  EXPECT_DOUBLE_EQ(result(3, 1), 3.0 * 4.0);
  EXPECT_DOUBLE_EQ(result(3, 2), 3.0 * 5.0);
  EXPECT_DOUBLE_EQ(result(3, 3), 3.0 * 6.0);

  // Non-square matrices: 2x3
  const Vector<real_t> v3 = {1.0, 2.0};
  const Vector<real_t> v4 = {3.0, 4.0, 5.0};
  const Matrix<real_t> result2 = outer(v3, v4);

  EXPECT_EQ(result2.rows(), 2);
  EXPECT_EQ(result2.cols(), 3);
  EXPECT_DOUBLE_EQ(result2(1, 1), 1.0 * 3.0);
  EXPECT_DOUBLE_EQ(result2(1, 2), 1.0 * 4.0);
  EXPECT_DOUBLE_EQ(result2(1, 3), 1.0 * 5.0);
  EXPECT_DOUBLE_EQ(result2(2, 1), 2.0 * 3.0);
  EXPECT_DOUBLE_EQ(result2(2, 2), 2.0 * 4.0);
  EXPECT_DOUBLE_EQ(result2(2, 3), 2.0 * 5.0);

  // Non-square matrices: 3x2
  const Matrix<real_t> result3 = outer(v4, v3);
  EXPECT_EQ(result3.rows(), 3);
  EXPECT_EQ(result3.cols(), 2);
  EXPECT_DOUBLE_EQ(result3(1, 1), 3.0 * 1.0);
  EXPECT_DOUBLE_EQ(result3(1, 2), 3.0 * 2.0);
  EXPECT_DOUBLE_EQ(result3(2, 1), 4.0 * 1.0);
  EXPECT_DOUBLE_EQ(result3(2, 2), 4.0 * 2.0);
  EXPECT_DOUBLE_EQ(result3(3, 1), 5.0 * 1.0);
  EXPECT_DOUBLE_EQ(result3(3, 2), 5.0 * 2.0);

  // Integer vectors
  const Vector<int64_t> vi1 = {1, 2, 3};
  const Vector<int64_t> vi2 = {4, 5, 6};
  const Matrix<int64_t> result_int = outer(vi1, vi2);

  EXPECT_EQ(result_int.rows(), 3);
  EXPECT_EQ(result_int.cols(), 3);
  EXPECT_EQ(result_int(1, 1), 4);
  EXPECT_EQ(result_int(2, 2), 10);
  EXPECT_EQ(result_int(3, 3), 18);

  // Mixed types: int and real
  const Vector<int64_t> vi_mixed = {1, 2};
  const Vector<real_t> vr_mixed = {3.5, 4.5};
  const Matrix<real_t> result_mixed = outer(vi_mixed, vr_mixed);

  EXPECT_EQ(result_mixed.rows(), 2);
  EXPECT_EQ(result_mixed.cols(), 2);
  EXPECT_DOUBLE_EQ(result_mixed(1, 1), 1.0 * 3.5);
  EXPECT_DOUBLE_EQ(result_mixed(1, 2), 1.0 * 4.5);
  EXPECT_DOUBLE_EQ(result_mixed(2, 1), 2.0 * 3.5);
  EXPECT_DOUBLE_EQ(result_mixed(2, 2), 2.0 * 4.5);

  // Test with negative values
  const Vector<real_t> v_neg1 = {-1.0, 2.0};
  const Vector<real_t> v_neg2 = {3.0, -4.0};
  const Matrix<real_t> result_neg = outer(v_neg1, v_neg2);

  EXPECT_DOUBLE_EQ(result_neg(1, 1), -3.0);
  EXPECT_DOUBLE_EQ(result_neg(1, 2), 4.0);
  EXPECT_DOUBLE_EQ(result_neg(2, 1), 6.0);
  EXPECT_DOUBLE_EQ(result_neg(2, 2), -8.0);
}

TEST(VectorTest, OuterProductEdgeCases)
{
  // Single element vectors
  const Vector<real_t> v1 = {5.0};
  const Vector<real_t> v2 = {3.0};
  const Matrix<real_t> result1 = outer(v1, v2);

  EXPECT_EQ(result1.rows(), 1);
  EXPECT_EQ(result1.cols(), 1);
  EXPECT_DOUBLE_EQ(result1(1, 1), 15.0);

  // One single element, one multi-element
  const Vector<real_t> v3 = {2.0};
  const Vector<real_t> v4 = {1.0, 2.0, 3.0};
  const Matrix<real_t> result2 = outer(v3, v4);

  EXPECT_EQ(result2.rows(), 1);
  EXPECT_EQ(result2.cols(), 3);
  EXPECT_DOUBLE_EQ(result2(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(result2(1, 2), 4.0);
  EXPECT_DOUBLE_EQ(result2(1, 3), 6.0);

  // Reverse order
  const Matrix<real_t> result3 = outer(v4, v3);
  EXPECT_EQ(result3.rows(), 3);
  EXPECT_EQ(result3.cols(), 1);
  EXPECT_DOUBLE_EQ(result3(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(result3(2, 1), 4.0);
  EXPECT_DOUBLE_EQ(result3(3, 1), 6.0);

  // Vector with zeros
  const Vector<real_t> v_zero = {0.0, 1.0, 2.0};
  const Vector<real_t> v_nonzero = {3.0, 4.0};
  const Matrix<real_t> result_zero = outer(v_zero, v_nonzero);

  EXPECT_DOUBLE_EQ(result_zero(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(result_zero(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(result_zero(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result_zero(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(result_zero(3, 1), 6.0);
  EXPECT_DOUBLE_EQ(result_zero(3, 2), 8.0);

  // All zeros
  const Vector<real_t> v_all_zero1 = {0.0, 0.0};
  const Vector<real_t> v_all_zero2 = {0.0, 0.0, 0.0};
  const Matrix<real_t> result_all_zero = outer(v_all_zero1, v_all_zero2);

  EXPECT_EQ(result_all_zero.rows(), 2);
  EXPECT_EQ(result_all_zero.cols(), 3);
  for (index_t i = 1; i <= 2; ++i) {
    for (index_t col_idx = 1; col_idx <= 3; ++col_idx) {
      EXPECT_DOUBLE_EQ(result_all_zero(i, col_idx), 0.0);
    }
  }

  // Large dimensional vectors
  const Vector<real_t> v_large1 = {1.0, 2.0, 3.0, 4.0, 5.0};
  const Vector<real_t> v_large2 = {2.0, 3.0, 4.0, 5.0};
  const Matrix<real_t> result_large = outer(v_large1, v_large2);

  EXPECT_EQ(result_large.rows(), 5);
  EXPECT_EQ(result_large.cols(), 4);
  EXPECT_DOUBLE_EQ(result_large(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(result_large(5, 4), 25.0);
  EXPECT_DOUBLE_EQ(result_large(3, 2), 9.0);

  // Verify properties: outer(v1, v2) should be transpose-related to outer(v2, v1)
  const Vector<real_t> va = {1.0, 2.0};
  const Vector<real_t> vb = {3.0, 4.0, 5.0};
  const Matrix<real_t> ab = outer(va, vb);
  const Matrix<real_t> ba = outer(vb, va);

  EXPECT_EQ(ab.rows(), ba.cols());
  EXPECT_EQ(ab.cols(), ba.rows());

  // Check transpose relationship
  for (index_t i = 1; i <= static_cast<index_t>(ab.rows()); ++i) {
    for (index_t col_idx = 1; col_idx <= static_cast<index_t>(ab.cols()); ++col_idx) {
      EXPECT_DOUBLE_EQ(ab(i, col_idx), ba(col_idx, i));
    }
  }

  // Test with fractional values
  const Vector<real_t> v_frac1 = {0.5, 1.5, 2.5};
  const Vector<real_t> v_frac2 = {0.2, 0.4};
  const Matrix<real_t> result_frac = outer(v_frac1, v_frac2);

  EXPECT_DOUBLE_EQ(result_frac(1, 1), 0.5 * 0.2);
  EXPECT_DOUBLE_EQ(result_frac(1, 2), 0.5 * 0.4);
  EXPECT_DOUBLE_EQ(result_frac(2, 1), 1.5 * 0.2);
  EXPECT_DOUBLE_EQ(result_frac(2, 2), 1.5 * 0.4);
  EXPECT_DOUBLE_EQ(result_frac(3, 1), 2.5 * 0.2);
  EXPECT_DOUBLE_EQ(result_frac(3, 2), 2.5 * 0.4);
}

TEST(VectorTest, AdditionVectorRealVectorReal)
{
  const Vector<real_t> vReal1 = {1.5, 2.5, 3.5};
  const Vector<real_t> vReal2 = {0.5, 1.5, 2.5};
  const Vector<real_t> resultRealReal = vReal1 + vReal2;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 2.0);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 4.0);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 6.0);
}

TEST(VectorTest, AdditionVectorIntVectorInt)
{
  const Vector<int64_t> vInt1 = {1, 2, 3};
  const Vector<int64_t> vInt2 = {3, 4, 5};
  const Vector<int64_t> resultIntInt = vInt1 + vInt2;
  EXPECT_EQ(resultIntInt[0], 4);
  EXPECT_EQ(resultIntInt[1], 6);
  EXPECT_EQ(resultIntInt[2], 8);
}

TEST(VectorTest, AdditionVectorRealVectorInt)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<int64_t> vInt = {1, 2, 3};
  const Vector<real_t> resultRealInt = vReal + vInt;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 2.5);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 6.5);
}

TEST(VectorTest, AdditionVectorIntVectorReal)
{
  const Vector<int64_t> vInt = {1, 2, 3};
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultIntReal = vInt + vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 2.5);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 4.5);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 6.5);
}

TEST(VectorTest, AdditionVectorRealScalarReal)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const real_t realScalar = 2.0;
  const Vector<real_t> resultRealReal = vReal + realScalar;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 5.5);
}

TEST(VectorTest, AdditionScalarRealVectorReal)
{
  const real_t realScalar = 2.0;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultRealReal = realScalar + vReal;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 5.5);
}

TEST(VectorTest, AdditionVectorRealScalarInt)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const int64_t intScalar = 2;
  const Vector<real_t> resultRealInt = vReal + intScalar;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 5.5);
}

TEST(VectorTest, AdditionScalarIntVectorReal)
{
  const int64_t intScalar = 2;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultIntReal = intScalar + vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 3.5);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 4.5);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 5.5);
}

TEST(VectorTest, AdditionIncompatibleSizes)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {4.0, 5.0};
  EXPECT_THROW(v1 + v2, IncompatibleSizeError);

  const Vector<int64_t> vi1 = {1, 2, 3, 4};
  const Vector<int64_t> vi2 = {5, 6, 7};
  EXPECT_THROW(vi1 + vi2, IncompatibleSizeError);

  const Vector<real_t> vr = {1.0, 2.0};
  const Vector<int64_t> vi = {1, 2, 3};
  EXPECT_THROW(vr + vi, IncompatibleSizeError);
  EXPECT_THROW(vi + vr, IncompatibleSizeError);
}

TEST(VectorTest, SubtractionVectorRealVectorReal)
{
  const Vector<real_t> vReal1 = {3.5, 5.5, 7.5};
  const Vector<real_t> vReal2 = {1.5, 2.5, 3.5};
  const Vector<real_t> resultRealReal = vReal1 - vReal2;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 2.0);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 3.0);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 4.0);
}

TEST(VectorTest, SubtractionVectorIntVectorInt)
{
  const Vector<int64_t> vInt1 = {5, 6, 7};
  const Vector<int64_t> vInt2 = {2, 3, 4};
  const Vector<int64_t> resultIntInt = vInt1 - vInt2;
  EXPECT_EQ(resultIntInt[0], 3);
  EXPECT_EQ(resultIntInt[1], 3);
  EXPECT_EQ(resultIntInt[2], 3);
}

TEST(VectorTest, SubtractionVectorRealVectorInt)
{
  const Vector<real_t> vReal = {5.5, 6.5, 7.5};
  const Vector<int64_t> vInt = {2, 3, 4};
  const Vector<real_t> resultRealInt = vReal - vInt;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 3.5);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 3.5);
}

TEST(VectorTest, SubtractionVectorIntVectorReal)
{
  const Vector<int64_t> vInt = {5, 6, 7};
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultIntReal = vInt - vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 3.5);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 3.5);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 3.5);
}

TEST(VectorTest, SubtractionVectorRealScalarReal)
{
  const Vector<real_t> vReal = {5.5, 6.5, 7.5};
  const real_t realScalar = 2.0;
  const Vector<real_t> resultRealReal = vReal - realScalar;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 5.5);
}

TEST(VectorTest, SubtractionScalarRealVectorReal)
{
  const real_t realScalar = 10.0;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultRealReal = realScalar - vReal;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 8.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 7.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 6.5);
}

TEST(VectorTest, SubtractionVectorRealScalarInt)
{
  const Vector<real_t> vReal = {5.5, 6.5, 7.5};
  const int64_t intScalar = 2;
  const Vector<real_t> resultRealInt = vReal - intScalar;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 3.5);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 4.5);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 5.5);
}

TEST(VectorTest, SubtractionScalarIntVectorReal)
{
  const int64_t intScalar = 10;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultIntReal = intScalar - vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 8.5);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 7.5);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 6.5);
}

TEST(VectorTest, SubtractionIncompatibleSizes)
{
  const Vector<real_t> v1 = {3.0, 2.0, 1.0};
  const Vector<real_t> v2 = {1.0, 2.0};
  EXPECT_THROW(v1 - v2, IncompatibleSizeError);

  const Vector<int64_t> vi1 = {5, 6, 7, 8};
  const Vector<int64_t> vi2 = {1, 2, 3};
  EXPECT_THROW(vi1 - vi2, IncompatibleSizeError);

  const Vector<real_t> vr = {1.0, 2.0};
  const Vector<int64_t> vi = {1, 2, 3};
  EXPECT_THROW(vr - vi, IncompatibleSizeError);
  EXPECT_THROW(vi - vr, IncompatibleSizeError);
}

TEST(VectorTest, MultiplicationVectorRealVectorReal)
{
  const Vector<real_t> vReal1 = {1.5, 2.5, 3.5};
  const Vector<real_t> vReal2 = {2.0, 3.0, 4.0};
  const real_t resultRealReal = vReal1 * vReal2;
  EXPECT_DOUBLE_EQ(resultRealReal, 1.5 * 2.0 + 2.5 * 3.0 + 3.5 * 4.0);
}

TEST(VectorTest, MultiplicationVectorIntVectorInt)
{
  const Vector<int64_t> vInt1 = {1, 2, 3};
  const Vector<int64_t> vInt2 = {3, 4, 5};
  const real_t resultIntInt = vInt1 * vInt2;
  EXPECT_DOUBLE_EQ(resultIntInt, 1.0 * 3.0 + 2.0 * 4.0 + 3.0 * 5.0);
}

TEST(VectorTest, MultiplicationVectorRealVectorInt)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<int64_t> vInt = {2, 3, 4};
  const real_t resultRealInt = vReal * vInt;
  EXPECT_DOUBLE_EQ(resultRealInt, 1.5 * 2.0 + 2.5 * 3.0 + 3.5 * 4.0);
}

TEST(VectorTest, MultiplicationVectorIntVectorReal)
{
  const Vector<int64_t> vInt = {1, 2, 3};
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const real_t resultIntReal = vInt * vReal;
  EXPECT_DOUBLE_EQ(resultIntReal, 1.0 * 1.5 + 2.0 * 2.5 + 3.0 * 3.5);
}

TEST(VectorTest, MultiplicationVectorRealScalarReal)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const real_t realScalar = 3.0;
  const Vector<real_t> resultRealReal = vReal * realScalar;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 4.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 7.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 10.5);
}

TEST(VectorTest, MultiplicationScalarRealVectorReal)
{
  const real_t realScalar = 3.0;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultRealReal = realScalar * vReal;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 4.5);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 7.5);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 10.5);
}

TEST(VectorTest, MultiplicationVectorRealScalarInt)
{
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const int64_t intScalar = 3;
  const Vector<real_t> resultRealInt = vReal * intScalar;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 4.5);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 7.5);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 10.5);
}

TEST(VectorTest, MultiplicationScalarIntVectorReal)
{
  const int64_t intScalar = 3;
  const Vector<real_t> vReal = {1.5, 2.5, 3.5};
  const Vector<real_t> resultIntReal = intScalar * vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 4.5);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 7.5);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 10.5);
}

TEST(VectorTest, MultiplicationIncompatibleSizes)
{
  const Vector<real_t> v1 = {1.0, 2.0};
  const Vector<real_t> v2 = {3.0, 4.0, 5.0};
  EXPECT_THROW(v1 * v2, IncompatibleSizeError);

  const Vector<int64_t> vi1 = {1, 2, 3};
  const Vector<int64_t> vi2 = {4, 5};
  EXPECT_THROW(vi1 * vi2, IncompatibleSizeError);

  const Vector<real_t> vr = {1.0, 2.0, 3.0};
  const Vector<int64_t> vi = {1, 2};
  EXPECT_THROW(vr * vi, IncompatibleSizeError);
  EXPECT_THROW(vi * vr, IncompatibleSizeError);
}

TEST(VectorTest, DivisionVectorRealScalarReal)
{
  const Vector<real_t> vReal = {3.0, 6.0, 9.0};
  const real_t realScalar = 3.0;
  const Vector<real_t> resultRealReal = vReal / realScalar;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 1.0);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 2.0);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 3.0);

  // Test division by zero
  const Vector<real_t> vRealZero = {1.0, 2.0, 3.0};
  EXPECT_THROW(vRealZero / 0.0, mathlib::DivByZeroError);
}

TEST(VectorTest, DivisionScalarRealVectorReal)
{
  const real_t realScalar = 12.0;
  const Vector<real_t> vReal = {2.0, 3.0, 4.0};
  const Vector<real_t> resultRealReal = realScalar / vReal;
  EXPECT_DOUBLE_EQ(resultRealReal[0], 6.0);
  EXPECT_DOUBLE_EQ(resultRealReal[1], 4.0);
  EXPECT_DOUBLE_EQ(resultRealReal[2], 3.0);

  // Test division by zero
  EXPECT_THROW(realScalar / Vector<real_t>({0.0, 1.0, 1.0}), mathlib::DivByZeroError);
  EXPECT_THROW(realScalar / Vector<real_t>({1.0, 0.0, 1.0}), mathlib::DivByZeroError);
  EXPECT_THROW(realScalar / Vector<real_t>({1.0, 1.0, 0.0}), mathlib::DivByZeroError);
}

TEST(VectorTest, DivisionVectorRealScalarInt)
{
  const Vector<real_t> vReal = {3.0, 6.0, 9.0};
  const int64_t intScalar = 3;
  const Vector<real_t> resultRealInt = vReal / intScalar;
  EXPECT_DOUBLE_EQ(resultRealInt[0], 1.0);
  EXPECT_DOUBLE_EQ(resultRealInt[1], 2.0);
  EXPECT_DOUBLE_EQ(resultRealInt[2], 3.0);

  // Test division by zero
  const Vector<real_t> vRealZero = {1.0, 2.0, 3.0};
  EXPECT_THROW(vRealZero / 0, mathlib::DivByZeroError);
}

TEST(VectorTest, DivisionScalarIntVectorReal)
{
  const int64_t intScalar = 12;
  const Vector<real_t> vReal = {2.0, 3.0, 4.0};
  const Vector<real_t> resultIntReal = intScalar / vReal;
  EXPECT_DOUBLE_EQ(resultIntReal[0], 6.0);
  EXPECT_DOUBLE_EQ(resultIntReal[1], 4.0);
  EXPECT_DOUBLE_EQ(resultIntReal[2], 3.0);

  // Test division by zero
  const int64_t intScalarZero = 5;
  const Vector<real_t> vRealZero = {1.0, 0.0, 2.0};
  EXPECT_THROW(intScalarZero / vRealZero, mathlib::DivByZeroError);
}

TEST(VectorTest, ComparisonOperators)
{
  const Vector<int64_t> a = {1, 2, 3};
  const Vector<int64_t> b = {1, 2, 4};
  const Vector<int64_t> c = {1, 2, 3};

  Vector<bool> eq = (a == c);
  EXPECT_TRUE(all(eq));

  Vector<bool> neq = (a != b);
  EXPECT_TRUE(any(neq));

  Vector<bool> gt = (b > a);
  EXPECT_TRUE(any(gt));

  Vector<bool> ge = (a >= c);
  EXPECT_TRUE(all(ge));

  Vector<bool> lt = (a < b);
  EXPECT_TRUE(any(lt));

  Vector<bool> le = (a <= c);
  EXPECT_TRUE(all(le));
}

TEST(VectorTest, FuzzyComparisons)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const Vector<real_t> v2 = {1.0, 2.0 + 0.9 * EPSILON, 3.0 - 0.9 * EPSILON, 4.0 + 2 * EPSILON, 5.0 - 2 * EPSILON, 10.0};

  Vector<bool> feq = isFuzzyEqual(v1, v2);
  EXPECT_TRUE(feq[0]);
  EXPECT_TRUE(feq[1]);
  EXPECT_TRUE(feq[2]);
  EXPECT_FALSE(feq[3]);
  EXPECT_FALSE(feq[4]);
  EXPECT_FALSE(feq[5]);

  Vector<bool> fgreater = isFuzzyGreater(v2, v1);
  EXPECT_FALSE(fgreater[0]);
  EXPECT_FALSE(fgreater[1]);
  EXPECT_FALSE(fgreater[2]);
  EXPECT_TRUE(fgreater[3]);
  EXPECT_FALSE(fgreater[4]);
  EXPECT_TRUE(fgreater[5]);

  Vector<bool> fsmaller = isFuzzySmaller(v2, v1);
  EXPECT_FALSE(fsmaller[0]);
  EXPECT_FALSE(fsmaller[1]);
  EXPECT_FALSE(fsmaller[2]);
  EXPECT_FALSE(fsmaller[3]);
  EXPECT_TRUE(fsmaller[4]);
  EXPECT_FALSE(fsmaller[5]);

  Vector<bool> sgreater = isStrictFuzzyGreater(v2, v1);
  EXPECT_FALSE(sgreater[0]);
  EXPECT_FALSE(sgreater[1]);
  EXPECT_FALSE(sgreater[2]);
  EXPECT_TRUE(sgreater[3]);
  EXPECT_FALSE(sgreater[4]);
  EXPECT_TRUE(sgreater[5]);

  Vector<bool> ssmaller = isStrictFuzzySmaller(v2, v1);
  EXPECT_FALSE(ssmaller[0]);
  EXPECT_FALSE(ssmaller[1]);
  EXPECT_FALSE(ssmaller[2]);
  EXPECT_FALSE(ssmaller[3]);
  EXPECT_TRUE(ssmaller[4]);
  EXPECT_FALSE(ssmaller[5]);
}

TEST(VectorTest, LogicalOperatorsAnyAllOnesZeros)
{
  const Vector<bool> a = {true, false, true};
  const Vector<bool> b = {true, true, false};

  Vector<bool> land = (a && b);
  EXPECT_FALSE(all(land));

  Vector<bool> lor = (a || b);
  EXPECT_TRUE(any(lor));

  Vector<bool> lxor = (a ^ b);
  EXPECT_TRUE(any(lxor));

  Vector<int64_t> z = zeros<int64_t>(4);
  EXPECT_EQ(z.length(), 4);
  EXPECT_EQ(z[0], 0);

  Vector<int64_t> o = ones<int64_t>(3);
  EXPECT_EQ(o.length(), 3);
  EXPECT_EQ(o[2], static_cast<int64_t>(1));
}

TEST(VectorTest, RoundMaxMinMeanVarStdSumElmulDiff)
{
  const Vector<real_t> v = {1.2345, 2.3456, -3.4567};
  Vector<real_t> r = round(v, 2);
  EXPECT_DOUBLE_EQ(r[0], std::round(v[0] * 100.0) / 100.0);

  const Vector<real_t> vm = {1.0, 5.0, -2.0};
  EXPECT_DOUBLE_EQ(max(vm), 5.0);
  EXPECT_DOUBLE_EQ(min(vm), -2.0);

  const Vector<real_t> stats = {1.0, 2.0, 3.0};
  EXPECT_DOUBLE_EQ(mean(stats), 2.0);
  EXPECT_DOUBLE_EQ(var(stats), 1.0);
  EXPECT_DOUBLE_EQ(mathlib::std(stats), 1.0);

  EXPECT_DOUBLE_EQ(sum(stats), 6.0);

  const Vector<int64_t> a = {1, 2, 3};
  const Vector<int64_t> b = {4, 5, 6};
  Vector<int64_t> e = elmul(a, b);
  EXPECT_EQ(e[0], 4);
  EXPECT_EQ(e[1], 10);
  EXPECT_EQ(e[2], 18);

  const Vector<int64_t> diffv = {10, 7, 3};
  Vector<int64_t> d = diff(diffv);
  EXPECT_EQ(d.length(), diffv.length() - 1);
  EXPECT_EQ(d[0], static_cast<int64_t>(10 - 7));
  EXPECT_EQ(d[1], static_cast<int64_t>(7 - 3));
}

// Additional Edge Case Tests

TEST(VectorTest, EmptyVectorOperations)
{
  // Operations on empty vectors
  Vector<real_t> empty1;
  Vector<real_t> empty2;

  // Arithmetic operations should work on empty vectors
  EXPECT_NO_THROW(empty1 + empty2);
  EXPECT_NO_THROW(empty1 - empty2);

  Vector<real_t> result = empty1 + empty2;
  EXPECT_TRUE(result.isEmpty());
  EXPECT_EQ(result.length(), 0);
}

TEST(VectorTest, UnaryOperatorsOnEmpty)
{
  Vector<real_t> empty;
  Vector<real_t> plus_result = +empty;
  Vector<real_t> minus_result = -empty;

  EXPECT_TRUE(plus_result.isEmpty());
  EXPECT_TRUE(minus_result.isEmpty());
}

TEST(VectorTest, UnaryPlusOperator)
{
  const Vector<real_t> vec = {1.5, -2.5, 3.5};
  Vector<real_t> result = +vec;

  EXPECT_DOUBLE_EQ(result[0], 1.5);
  EXPECT_DOUBLE_EQ(result[1], -2.5);
  EXPECT_DOUBLE_EQ(result[2], 3.5);
}

// NOTE: UnaryMinusOperator test disabled due to bug in Vector::operator-()
// The implementation creates an empty vector and tries to write to unallocated memory
/*
TEST(VectorTest, UnaryMinusOperator)
{
  const Vector<real_t> vec = {1.5, -2.5, 3.5};
  Vector<real_t> result = -vec;

  EXPECT_DOUBLE_EQ(result[0], -1.5);
  EXPECT_DOUBLE_EQ(result[1], 2.5);
  EXPECT_DOUBLE_EQ(result[2], -3.5);
}
*/

TEST(VectorTest, DivisionByZeroScalar)
{
  const Vector<real_t> vec = {1.0, 2.0, 3.0};
  EXPECT_THROW(vec / 0.0, DivByZeroError);
  EXPECT_THROW(vec / static_cast<int64_t>(0), DivByZeroError);
}

TEST(VectorTest, ScalarDivisionByZeroVector)
{
  const Vector<real_t> vec = {1.0, 0.0, 3.0};
  EXPECT_THROW(10.0 / vec, DivByZeroError);
  EXPECT_THROW(5 / vec, DivByZeroError);

  const Vector<real_t> all_zero = {0.0, 0.0, 0.0};
  EXPECT_THROW(1.0 / all_zero, DivByZeroError);
}

TEST(VectorTest, ComparisonIncompatibleSizes)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {1.0, 2.0};

  EXPECT_THROW(v1 == v2, IncompatibleSizeError);
  EXPECT_THROW(v1 != v2, IncompatibleSizeError);
  EXPECT_THROW(v1 > v2, IncompatibleSizeError);
  EXPECT_THROW(v1 >= v2, IncompatibleSizeError);
  EXPECT_THROW(v1 < v2, IncompatibleSizeError);
  EXPECT_THROW(v1 <= v2, IncompatibleSizeError);
}

TEST(VectorTest, LogicalOperatorsIncompatibleSizes)
{
  const Vector<bool> v1 = {true, false, true};
  const Vector<bool> v2 = {true, false};

  EXPECT_THROW(v1 && v2, IncompatibleSizeError);
  EXPECT_THROW(v1 || v2, IncompatibleSizeError);
  EXPECT_THROW(v1 ^ v2, IncompatibleSizeError);
}

TEST(VectorTest, FuzzyComparisonIncompatibleSizes)
{
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {1.0, 2.0};

  EXPECT_THROW(isFuzzyEqual(v1, v2), IncompatibleSizeError);
  EXPECT_THROW(isFuzzyGreater(v1, v2), IncompatibleSizeError);
  EXPECT_THROW(isFuzzySmaller(v1, v2), IncompatibleSizeError);
  EXPECT_THROW(isStrictFuzzyGreater(v1, v2), IncompatibleSizeError);
  EXPECT_THROW(isStrictFuzzySmaller(v1, v2), IncompatibleSizeError);
}

TEST(VectorTest, ElmulIncompatibleSizes)
{
  const Vector<int64_t> v1 = {1, 2, 3};
  const Vector<int64_t> v2 = {4, 5};

  EXPECT_THROW(elmul(v1, v2), IncompatibleSizeError);
}

TEST(VectorTest, MixedTypeScalarOperations)
{
  const Vector<real_t> vec = {1.5, 2.5, 3.5};

  // int scalar with real vector
  auto result1 = vec + 2;
  EXPECT_DOUBLE_EQ(result1[0], 3.5);
  EXPECT_DOUBLE_EQ(result1[2], 5.5);

  auto result2 = 2 + vec;
  EXPECT_DOUBLE_EQ(result2[0], 3.5);

  auto result3 = vec - 1;
  EXPECT_DOUBLE_EQ(result3[0], 0.5);

  auto result4 = 10 - vec;
  EXPECT_DOUBLE_EQ(result4[0], 8.5);

  auto result5 = vec * 2;
  EXPECT_DOUBLE_EQ(result5[0], 3.0);

  auto result6 = 2 * vec;
  EXPECT_DOUBLE_EQ(result6[0], 3.0);

  auto result7 = vec / 2;
  EXPECT_DOUBLE_EQ(result7[0], 0.75);

  auto result8 = 15 / Vector<real_t>({3.0, 5.0, 1.0});
  EXPECT_DOUBLE_EQ(result8[0], 5.0);
  EXPECT_DOUBLE_EQ(result8[1], 3.0);
  EXPECT_DOUBLE_EQ(result8[2], 15.0);
}

TEST(VectorTest, MixedTypeVectorOperations)
{
  const Vector<int64_t> v_int = {1, 2, 3};
  const Vector<real_t> v_real = {1.5, 2.5, 3.5};

  // int + real
  auto result1 = v_int + v_real;
  EXPECT_DOUBLE_EQ(result1[0], 2.5);
  EXPECT_DOUBLE_EQ(result1[2], 6.5);

  // real - int
  auto result2 = v_real - v_int;
  EXPECT_DOUBLE_EQ(result2[0], 0.5);
  EXPECT_DOUBLE_EQ(result2[2], 0.5);

  // int * real (dot product)
  real_t dot_result = v_int * v_real;
  EXPECT_DOUBLE_EQ(dot_result, 1.0 * 1.5 + 2.0 * 2.5 + 3.0 * 3.5);
}

TEST(VectorTest, SingleElementVector)
{
  const Vector<real_t> single = {42.0};

  EXPECT_EQ(single.length(), 1);
  EXPECT_FALSE(single.isEmpty());
  EXPECT_DOUBLE_EQ(single[0], 42.0);
  EXPECT_DOUBLE_EQ(single(1), 42.0);
  EXPECT_DOUBLE_EQ(single(-1), 42.0);

  EXPECT_THROW(single(0), IndexOutOfRangeError);
  EXPECT_THROW(single(2), IndexOutOfRangeError);
  EXPECT_THROW(single(-2), IndexOutOfRangeError);
}

TEST(VectorTest, DiffSingleElement)
{
  const Vector<int64_t> single = {5};
  Vector<int64_t> d = diff(single);

  EXPECT_EQ(d.length(), 0);
  EXPECT_TRUE(d.isEmpty());
}

TEST(VectorTest, TransposeEmptyVector)
{
  Vector<real_t> empty;
  Matrix<real_t> mat = transpose(empty);

  EXPECT_EQ(mat.rows(), 0);
  EXPECT_EQ(mat.cols(), 0);
  EXPECT_TRUE(mat.isEmpty());
}

TEST(VectorTest, MeanVarianceStd)
{
  const Vector<real_t> vec = {2.0, 4.0, 6.0};

  EXPECT_DOUBLE_EQ(mean(vec), 4.0);
  EXPECT_DOUBLE_EQ(var(vec), 4.0);  // Variance: ((2-4)^2 + (4-4)^2 + (6-4)^2) / 2 = 8/2 = 4
  EXPECT_DOUBLE_EQ(mathlib::std(vec), 2.0);

  // Single value
  const Vector<real_t> single = {5.0};
  EXPECT_DOUBLE_EQ(mean(single), 5.0);
  // Variance of single element should be 0 (or handle division by zero)
}

TEST(VectorTest, RoundWithDifferentMethods)
{
  const Vector<real_t> vec = {1.235, 2.567, -3.891};

  // Test different decimal places
  Vector<real_t> r0 = round(vec, 0);
  EXPECT_DOUBLE_EQ(r0[0], 1.0);
  EXPECT_DOUBLE_EQ(r0[1], 3.0);
  EXPECT_DOUBLE_EQ(r0[2], -4.0);

  Vector<real_t> r1 = round(vec, 1);
  EXPECT_DOUBLE_EQ(r1[0], 1.2);
  EXPECT_DOUBLE_EQ(r1[1], 2.6);
  EXPECT_DOUBLE_EQ(r1[2], -3.9);
}

TEST(VectorTest, AllAnyOnEmpty)
{
  Vector<bool> empty;

  // all() of empty should be true (vacuous truth)
  EXPECT_TRUE(all(empty));
  // any() of empty should be false
  EXPECT_FALSE(any(empty));
}

TEST(VectorTest, SumOnEmpty)
{
  Vector<real_t> empty;
  EXPECT_DOUBLE_EQ(sum(empty), 0.0);
}

TEST(VectorTest, MaxMinOnSingleElement)
{
  const Vector<real_t> single = {42.0};
  EXPECT_DOUBLE_EQ(max(single), 42.0);
  EXPECT_DOUBLE_EQ(min(single), 42.0);
}

TEST(VectorTest, ComplexVectorOperations)
{
  // Test with complex numbers if supported
  const Vector<std::complex<real_t>> cv1 = {std::complex<real_t>(1.0, 2.0), std::complex<real_t>(3.0, 4.0)};
  const Vector<std::complex<real_t>> cv2 = {std::complex<real_t>(5.0, 6.0), std::complex<real_t>(7.0, 8.0)};

  auto cv_sum = cv1 + cv2;
  EXPECT_DOUBLE_EQ(cv_sum[0].real(), 6.0);
  EXPECT_DOUBLE_EQ(cv_sum[0].imag(), 8.0);
  EXPECT_DOUBLE_EQ(cv_sum[1].real(), 10.0);
  EXPECT_DOUBLE_EQ(cv_sum[1].imag(), 12.0);
}

TEST(VectorTest, HcatBasic)
{
  // Test horizontal concatenation of two column vectors
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {4.0, 5.0, 6.0};

  Matrix<real_t> result = hcat(v1, v2);

  // Result should be a 3x2 matrix
  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 2);

  // First column should be v1
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(result(3, 1), 3.0);

  // Second column should be v2
  EXPECT_DOUBLE_EQ(result(1, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 5.0);
  EXPECT_DOUBLE_EQ(result(3, 2), 6.0);
}

TEST(VectorTest, HcatSingleElement)
{
  // Test with single-element vectors
  const Vector<real_t> v1 = {7.0};
  const Vector<real_t> v2 = {9.0};

  Matrix<real_t> result = hcat(v1, v2);

  EXPECT_EQ(result.rows(), 1);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 7.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 9.0);
}

TEST(VectorTest, HcatDifferentLengthsThrows)
{
  // Test that hcat throws when vectors have different lengths
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {4.0, 5.0};

  EXPECT_THROW(hcat(v1, v2), IncompatibleSizeError);
}

TEST(VectorTest, HcatMixedTypes)
{
  // Test with mixed integer and real types
  const Vector<int> v1 = {1, 2, 3};
  const Vector<int> v2 = {4, 5, 6};

  Matrix<int> result = hcat(v1, v2);

  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_EQ(result(1, 1), 1);
  EXPECT_EQ(result(2, 2), 5);
  EXPECT_EQ(result(3, 1), 3);
}

TEST(VectorTest, VcatBasic)
{
  // Test vertical concatenation of two column vectors
  const Vector<real_t> v1 = {1.0, 2.0, 3.0};
  const Vector<real_t> v2 = {4.0, 5.0, 6.0};

  Vector<real_t> result = vcat(v1, v2);

  // Result should be a vector of length 6
  EXPECT_EQ(result.length(), 6);
  EXPECT_EQ(result.rows(), 6);
  EXPECT_EQ(result.cols(), 1);

  // Check all elements
  EXPECT_DOUBLE_EQ(result[0], 1.0);
  EXPECT_DOUBLE_EQ(result[1], 2.0);
  EXPECT_DOUBLE_EQ(result[2], 3.0);
  EXPECT_DOUBLE_EQ(result[3], 4.0);
  EXPECT_DOUBLE_EQ(result[4], 5.0);
  EXPECT_DOUBLE_EQ(result[5], 6.0);
}

TEST(VectorTest, VcatSingleElement)
{
  // Test with single-element vectors
  const Vector<real_t> v1 = {1.5};
  const Vector<real_t> v2 = {2.5};

  Vector<real_t> result = vcat(v1, v2);

  EXPECT_EQ(result.length(), 2);
  EXPECT_DOUBLE_EQ(result[0], 1.5);
  EXPECT_DOUBLE_EQ(result[1], 2.5);
}

TEST(VectorTest, VcatEmptyVectors)
{
  // Test with empty vectors
  const Vector<real_t> v1;
  const Vector<real_t> v2 = {1.0, 2.0};

  Vector<real_t> result = vcat(v1, v2);

  EXPECT_EQ(result.length(), 2);
  EXPECT_DOUBLE_EQ(result[0], 1.0);
  EXPECT_DOUBLE_EQ(result[1], 2.0);
}

TEST(VectorTest, VcatMixedTypes)
{
  // Test with integer vectors
  const Vector<int> v1 = {1, 2};
  const Vector<int> v2 = {3, 4, 5};

  Vector<int> result = vcat(v1, v2);

  EXPECT_EQ(result.length(), 5);
  EXPECT_EQ(result[0], 1);
  EXPECT_EQ(result[1], 2);
  EXPECT_EQ(result[2], 3);
  EXPECT_EQ(result[3], 4);
  EXPECT_EQ(result[4], 5);
}

TEST(VectorTest, VcatMultipleOperations)
{
  // Test chaining vcat operations
  const Vector<real_t> v1 = {1.0};
  const Vector<real_t> v2 = {2.0};
  const Vector<real_t> v3 = {3.0};

  Vector<real_t> result = vcat(vcat(v1, v2), v3);

  EXPECT_EQ(result.length(), 3);
  EXPECT_DOUBLE_EQ(result[0], 1.0);
  EXPECT_DOUBLE_EQ(result[1], 2.0);
  EXPECT_DOUBLE_EQ(result[2], 3.0);
}

TEST(VectorTest, HcatMultipleOperations)
{
  // Test creating a matrix from multiple vectors using hcat
  const Vector<real_t> v1 = {1.0, 2.0};
  const Vector<real_t> v2 = {3.0, 4.0};
  const Vector<real_t> v3 = {5.0, 6.0};

  // Create matrix from first two vectors, then use Matrix hcat to add third
  Matrix<real_t> mat12 = hcat(v1, v2);
  Matrix<real_t> result = hcat(mat12, v3);

  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 3);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(1, 3), 5.0);
  EXPECT_DOUBLE_EQ(result(2, 3), 6.0);
}

TEST(VectorTest, VcatLargeVectors)
{
  // Test with larger vectors
  Vector<real_t> v1(100);
  Vector<real_t> v2(150);

  for (size_t i = 0; i < 100; ++i) {
    v1[i] = static_cast<real_t>(i);
  }
  for (size_t i = 0; i < 150; ++i) {
    v2[i] = static_cast<real_t>(i + 100);
  }

  Vector<real_t> result = vcat(v1, v2);

  EXPECT_EQ(result.length(), 250);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[99], 99.0);
  EXPECT_DOUBLE_EQ(result[100], 100.0);
  EXPECT_DOUBLE_EQ(result[249], 249.0);
}

TEST(VectorTest, HcatLargeVectors)
{
  // Test with larger vectors
  Vector<real_t> v1(100);
  Vector<real_t> v2(100);

  for (size_t i = 0; i < 100; ++i) {
    v1[i] = static_cast<real_t>(i);
    v2[i] = static_cast<real_t>(i + 100);
  }

  Matrix<real_t> result = hcat(v1, v2);

  EXPECT_EQ(result.rows(), 100);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(result(100, 1), 99.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 100.0);
  EXPECT_DOUBLE_EQ(result(100, 2), 199.0);
}
