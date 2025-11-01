/*****************************************************************************************
 * @file: Matrix_test.cpp
 *
 * @brief: Unit tests for Matrix class
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "Matrix.hpp"

#include <gtest/gtest.h>

#include "config.hpp"

using namespace mathlib;

TEST(MatrixTest, DefaultConstructor)
{
  Matrix<real_t> mat;
  EXPECT_EQ(mat.rows(), 0);
  EXPECT_EQ(mat.cols(), 0);
  EXPECT_EQ(mat.length(), 0);
  EXPECT_TRUE(mat.isEmpty());
  EXPECT_THROW(mat[0], IndexOutOfRangeError);
  EXPECT_THROW(mat(1), IndexOutOfRangeError);
  EXPECT_THROW(mat(-1), IndexOutOfRangeError);
}

TEST(MatrixTest, SizeConstructor)
{
  Matrix<real_t> mat(3, 4);
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 4);
  EXPECT_EQ(mat.length(), 12);
  EXPECT_FALSE(mat.isEmpty());
}

TEST(MatrixTest, InitializerListOfElementsConstructor)
{
  Matrix<real_t> mat({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
  EXPECT_EQ(mat.rows(), 2);
  EXPECT_EQ(mat.cols(), 3);
  EXPECT_DOUBLE_EQ(mat[0], 1.0);
  EXPECT_DOUBLE_EQ(mat[1], 4.0);
  EXPECT_DOUBLE_EQ(mat[2], 2.0);
  EXPECT_DOUBLE_EQ(mat[3], 5.0);
  EXPECT_DOUBLE_EQ(mat[4], 3.0);
  EXPECT_DOUBLE_EQ(mat[5], 6.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 3), 3.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 4.0);
  EXPECT_DOUBLE_EQ(mat(2, 2), 5.0);
  EXPECT_DOUBLE_EQ(mat(2, 3), 6.0);
}

TEST(MatrixTest, InitializerListOfVectorsConstructor)
{
  Vector<real_t> vec1({1.0, 2.0, 3.0});
  Vector<real_t> vec2({4.0, 5.0, 6.0});
  Matrix<real_t> mat({vec1, vec2});
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 2);
  EXPECT_DOUBLE_EQ(mat[0], 1.0);
  EXPECT_DOUBLE_EQ(mat[1], 2.0);
  EXPECT_DOUBLE_EQ(mat[2], 3.0);
  EXPECT_DOUBLE_EQ(mat[3], 4.0);
  EXPECT_DOUBLE_EQ(mat[4], 5.0);
  EXPECT_DOUBLE_EQ(mat[5], 6.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 4.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(2, 2), 5.0);
  EXPECT_DOUBLE_EQ(mat(3, 1), 3.0);
  EXPECT_DOUBLE_EQ(mat(3, 2), 6.0);
}

TEST(MatrixTest, VectorConstructor)
{
  Vector<real_t> vec = {1.0, 2.0, 3.0};

  Matrix<real_t> mat(vec);
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 1);
  EXPECT_DOUBLE_EQ(mat[0], 1.0);
  EXPECT_DOUBLE_EQ(mat[1], 2.0);
  EXPECT_DOUBLE_EQ(mat[2], 3.0);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(3, 1), 3.0);
}

TEST(MatrixTest, CopyConstructor)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2(mat1);
  EXPECT_EQ(mat2.rows(), 2);
  EXPECT_EQ(mat2.cols(), 2);
  for (size_t i = 0; i < mat1.length(); ++i) {
    EXPECT_DOUBLE_EQ(mat2[i], mat1[i]);
  }
}

TEST(MatrixTest, MoveConstructor)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2(std::move(mat1));
  EXPECT_EQ(mat2.rows(), 2);
  EXPECT_EQ(mat2.cols(), 2);
  EXPECT_EQ(mat1.rows(), 0);
  EXPECT_EQ(mat1.cols(), 0);
}

TEST(MatrixTest, CopyAssignment)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2;
  mat2 = mat1;
  EXPECT_EQ(mat2.rows(), 2);
  EXPECT_EQ(mat2.cols(), 2);
  for (size_t i = 0; i < mat1.length(); ++i) {
    EXPECT_DOUBLE_EQ(mat2[i], mat1[i]);
  }
}

TEST(MatrixTest, MoveAssignment)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2;
  mat2 = std::move(mat1);
  EXPECT_EQ(mat2.rows(), 2);
  EXPECT_EQ(mat2.cols(), 2);
  EXPECT_EQ(mat1.rows(), 0);
  EXPECT_EQ(mat1.cols(), 0);
}

TEST(MatrixTest, LinearIndexing0Based)
{
  Matrix<real_t> mat({{1.0, 3.0, 5.0}, {2.0, 4.0, 6.0}});
  EXPECT_DOUBLE_EQ(mat[0], 1.0);
  EXPECT_DOUBLE_EQ(mat[1], 2.0);
  EXPECT_DOUBLE_EQ(mat[2], 3.0);
  EXPECT_DOUBLE_EQ(mat[3], 4.0);
  EXPECT_DOUBLE_EQ(mat[4], 5.0);
  EXPECT_DOUBLE_EQ(mat[5], 6.0);
  EXPECT_THROW(mat[6], IndexOutOfRangeError);
}

TEST(MatrixTest, LinearIndexing1Based)
{
  Matrix<real_t> mat({{1.0, 3.0, 5.0}, {2.0, 4.0, 6.0}});
  EXPECT_THROW(mat(0), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(mat(1), 1.0);
  EXPECT_DOUBLE_EQ(mat(2), 2.0);
  EXPECT_DOUBLE_EQ(mat(3), 3.0);
  EXPECT_DOUBLE_EQ(mat(4), 4.0);
  EXPECT_DOUBLE_EQ(mat(5), 5.0);
  EXPECT_DOUBLE_EQ(mat(6), 6.0);
  EXPECT_THROW(mat(7), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(mat(-6), 1.0);
  EXPECT_DOUBLE_EQ(mat(-5), 2.0);
  EXPECT_DOUBLE_EQ(mat(-4), 3.0);
  EXPECT_DOUBLE_EQ(mat(-3), 4.0);
  EXPECT_DOUBLE_EQ(mat(-2), 5.0);
  EXPECT_DOUBLE_EQ(mat(-1), 6.0);
  EXPECT_THROW(mat(-7), IndexOutOfRangeError);
}

TEST(MatrixTest, DoubleIndexing1Based)
{
  Matrix<real_t> mat({{1.0, 3.0, 5.0}, {2.0, 4.0, 6.0}});
  EXPECT_THROW(mat(0, 1), IndexOutOfRangeError);
  EXPECT_THROW(mat(1, 0), IndexOutOfRangeError);
  EXPECT_THROW(mat(0, 0), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 2.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 3.0);
  EXPECT_DOUBLE_EQ(mat(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(mat(1, 3), 5.0);
  EXPECT_DOUBLE_EQ(mat(2, 3), 6.0);
  EXPECT_THROW(mat(3, 1), IndexOutOfRangeError);
  EXPECT_THROW(mat(1, 4), IndexOutOfRangeError);
  EXPECT_THROW(mat(3, 4), IndexOutOfRangeError);
  EXPECT_DOUBLE_EQ(mat(-1, -1), 6.0);
  EXPECT_DOUBLE_EQ(mat(-2, -1), 5.0);
  EXPECT_DOUBLE_EQ(mat(-1, -2), 4.0);
  EXPECT_DOUBLE_EQ(mat(-2, -2), 3.0);
  EXPECT_DOUBLE_EQ(mat(-1, -3), 2.0);
  EXPECT_DOUBLE_EQ(mat(-2, -3), 1.0);
  EXPECT_THROW(mat(-3, 1), IndexOutOfRangeError);
  EXPECT_THROW(mat(1, -4), IndexOutOfRangeError);
  EXPECT_THROW(mat(-3, -4), IndexOutOfRangeError);
}

TEST(MatrixTest, ZerosConstructor)
{
  Matrix<real_t> mat = zeros<real_t>(3, 2);
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 2);
  for (size_t i = 0; i < mat.length(); ++i) {
    EXPECT_DOUBLE_EQ(mat[i], 0.0);
  }
  for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      EXPECT_DOUBLE_EQ(mat(r, c), 0.0);
    }
  }
}

TEST(MatrixTest, OnesConstructor)
{
  Matrix<real_t> mat = ones<real_t>(2, 4);
  EXPECT_EQ(mat.rows(), 2);
  EXPECT_EQ(mat.cols(), 4);
  for (size_t i = 0; i < mat.length(); ++i) {
    EXPECT_DOUBLE_EQ(mat[i], 1.0);
  }
  for (index_t r = 1; r <= static_cast<index_t>(mat.rows()); ++r) {
    for (index_t c = 1; c <= static_cast<index_t>(mat.cols()); ++c) {
      EXPECT_DOUBLE_EQ(mat(r, c), 1.0);
    }
  }
}

TEST(MatrixTest, EyeConstructor)
{
  Matrix<real_t> mat = eye<real_t>(3);
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 3);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat(1, 3), 0.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(mat(2, 3), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 3), 1.0);
}

TEST(MatrixTest, DiagFromVector)
{
  Vector<real_t> vec = {1.0, 2.0, 3.0};
  Matrix<real_t> mat = diag(vec);
  EXPECT_EQ(mat.rows(), 3);
  EXPECT_EQ(mat.cols(), 3);
  EXPECT_DOUBLE_EQ(mat(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat(1, 3), 0.0);
  EXPECT_DOUBLE_EQ(mat(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat(2, 2), 2.0);
  EXPECT_DOUBLE_EQ(mat(2, 3), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 1), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 2), 0.0);
  EXPECT_DOUBLE_EQ(mat(3, 3), 3.0);
}

TEST(MatrixTest, DiagFromMatrix)
{
  Matrix<real_t> mat({{1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}});
  Vector<real_t> vec = diag(mat);
  EXPECT_EQ(vec.length(), 3);
  EXPECT_DOUBLE_EQ(vec[0], 1.0);
  EXPECT_DOUBLE_EQ(vec[1], 5.0);
  EXPECT_DOUBLE_EQ(vec[2], 9.0);
}

TEST(MatrixTest, AdditionMatrixMatrix)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0}, {7.0, 8.0}});
  Matrix<real_t> result = mat1 + mat2;
  EXPECT_DOUBLE_EQ(result(1, 1), 6.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 8.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 10.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 12.0);
}

TEST(MatrixTest, AdditionMatrixScalar)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> result = mat + 10.0;
  EXPECT_DOUBLE_EQ(result(1, 1), 11.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 12.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 13.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 14.0);
}

TEST(MatrixTest, AdditionScalarMatrix)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> result = 10.0 + mat;
  EXPECT_DOUBLE_EQ(result(1, 1), 11.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 12.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 13.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 14.0);
}

TEST(MatrixTest, SubtractionMatrixMatrix)
{
  Matrix<real_t> mat1({{5.0, 6.0}, {7.0, 8.0}});
  Matrix<real_t> mat2({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> result = mat1 - mat2;
  EXPECT_DOUBLE_EQ(result(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 4.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
}

TEST(MatrixTest, SubtractionMatrixScalar)
{
  Matrix<real_t> mat({{11.0, 12.0}, {13.0, 14.0}});
  Matrix<real_t> result = mat - 10.0;
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
}

TEST(MatrixTest, SubtractionScalarMatrix)
{
  Matrix<real_t> mat({{11.0, 12.0}, {13.0, 14.0}});
  Matrix<real_t> result = 10.0 - mat;
  EXPECT_DOUBLE_EQ(result(1, 1), -1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), -2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), -3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), -4.0);
}

TEST(MatrixTest, MultiplicationMatrixEye)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> I = eye(2);
  Matrix<real_t> result = mat * I;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
}

TEST(MatrixTest, MultiplicationEyeMatrix)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> I = eye(2);
  Matrix<real_t> result = I * mat;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
}

TEST(MatrixTest, MultiplicationMatrixMatrix)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0, 7.0}, {8.0, 9.0, 10.0}});
  Matrix<real_t> result = mat1 * mat2;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 3);
  EXPECT_DOUBLE_EQ(result(1, 1), 21.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 24.0);
  EXPECT_DOUBLE_EQ(result(1, 3), 27.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 47.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 54.0);
  EXPECT_DOUBLE_EQ(result(2, 3), 61.0);
}

TEST(MatrixTest, MultiplicationMatrixVector)
{
  Matrix<real_t> mat({{1.0, 3.0, 5.0}, {2.0, 4.0, 6.0}});
  Vector<real_t> vec = {7.0, 8.0, 9.0};
  Vector<real_t> result = mat * vec;
  EXPECT_EQ(result.length(), 2);
  EXPECT_DOUBLE_EQ(result(1), 76.0);
  EXPECT_DOUBLE_EQ(result(2), 100.0);
}

TEST(MatrixTest, MultiplicationVectorMatrix)
{
  Vector<real_t> vec = {1.0, 2.0, 3.0};
  Matrix<real_t> mat({{4.0, 5.0, 6.0}});
  Matrix<real_t> result = vec * mat;
  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 3);
  EXPECT_DOUBLE_EQ(result(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 5.0);
  EXPECT_DOUBLE_EQ(result(1, 3), 6.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 8.0);
  EXPECT_DOUBLE_EQ(result(3, 1), 12.0);
}

TEST(MatrixTest, MultiplicationMatrixScalar)
{
  Matrix<real_t> mat({{1.0, 3.0}, {2.0, 4.0}});
  Matrix<real_t> result = mat * 2.0;
  EXPECT_DOUBLE_EQ(result(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 6.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 4.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 8.0);
}

TEST(MatrixTest, DivisionMatrixScalar)
{
  Matrix<real_t> mat({{10.0, 30.0}, {20.0, 40.0}});
  Matrix<real_t> result = mat / 2.0;
  EXPECT_DOUBLE_EQ(result(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 15.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 10.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 20.0);
}

TEST(MatrixTest, UnaryPlus)
{
  Matrix<real_t> mat({{1.0, 3.0}, {2.0, 4.0}});
  Matrix<real_t> result = +mat;
  for (size_t i = 0; i < mat.length(); ++i) {
    EXPECT_DOUBLE_EQ(result[i], mat[i]);
  }
}

TEST(MatrixTest, UnaryMinus)
{
  Matrix<real_t> mat{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<real_t> result = -mat;
  for (size_t i = 0; i < mat.length(); ++i) {
    EXPECT_DOUBLE_EQ(result[i], -mat[i]);
  }
}

TEST(MatrixTest, ComparisonOperators)
{
  Matrix<int64_t> mat1({{1, 2}, {3, 4}});
  Matrix<int64_t> mat2({{1, 3}, {2, -4}});

  Matrix<bool> eq = (mat1 == mat2);
  EXPECT_TRUE(eq(1, 1));
  EXPECT_FALSE(eq(1, 2));
  EXPECT_FALSE(eq(2, 1));
  EXPECT_FALSE(eq(2, 2));

  Matrix<bool> gt = (mat1 > mat2);
  EXPECT_FALSE(gt(1, 1));
  EXPECT_FALSE(gt(1, 2));
  EXPECT_TRUE(gt(2, 1));
  EXPECT_TRUE(gt(2, 2));

  Matrix<bool> gteq = (mat1 >= mat2);
  EXPECT_TRUE(gteq(1, 1));
  EXPECT_FALSE(gteq(1, 2));
  EXPECT_TRUE(gteq(2, 1));
  EXPECT_TRUE(gteq(2, 2));

  Matrix<bool> lt = (mat1 < mat2);
  EXPECT_FALSE(lt(1, 1));
  EXPECT_TRUE(lt(1, 2));
  EXPECT_FALSE(lt(2, 1));
  EXPECT_FALSE(lt(2, 2));

  Matrix<bool> lteq = (mat1 <= mat2);
  EXPECT_TRUE(lteq(1, 1));
  EXPECT_TRUE(lteq(1, 2));
  EXPECT_FALSE(lteq(2, 1));
  EXPECT_FALSE(lteq(2, 2));
}

TEST(MatrixTest, LogicalOperators)
{
  Matrix<bool> mat1({{true, false}, {true, false}});
  Matrix<bool> mat2({{true, true}, {false, false}});

  Matrix<bool> andResult = (mat1 && mat2);
  EXPECT_TRUE(andResult(1, 1));
  EXPECT_FALSE(andResult(1, 2));
  EXPECT_FALSE(andResult(2, 1));
  EXPECT_FALSE(andResult(2, 2));

  Matrix<bool> orResult = (mat1 || mat2);
  EXPECT_TRUE(orResult(1, 1));
  EXPECT_TRUE(orResult(1, 2));
  EXPECT_TRUE(orResult(2, 1));
  EXPECT_FALSE(orResult(2, 2));

  Matrix<bool> exorResult = (mat1 ^ mat2);
  EXPECT_FALSE(exorResult(1, 1));
  EXPECT_TRUE(exorResult(1, 2));
  EXPECT_TRUE(exorResult(2, 1));
  EXPECT_FALSE(exorResult(2, 2));
}

TEST(MatrixTest, Transpose)
{
  Matrix<real_t> mat({{1.0, 3.0, 5.0}, {2.0, 4.0, 6.0}});
  EXPECT_EQ(mat.rows(), 2);
  EXPECT_EQ(mat.cols(), 3);
  Matrix<real_t> result = transpose(mat);
  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(3, 1), 5.0);
  EXPECT_DOUBLE_EQ(result(3, 2), 6.0);
}

TEST(MatrixTest, Trace)
{
  Matrix<real_t> mat({{1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}});
  real_t result = trace(mat);
  EXPECT_DOUBLE_EQ(result, 1.0 + 5.0 + 9.0);
}

TEST(MatrixTest, Abs)
{
  Matrix<real_t> mat({{-1.0, -3.0}, {2.0, 4.0}});
  Matrix<real_t> result = abs(mat);
  EXPECT_DOUBLE_EQ(result[0], 1.0);
  EXPECT_DOUBLE_EQ(result[1], 2.0);
  EXPECT_DOUBLE_EQ(result[2], 3.0);
  EXPECT_DOUBLE_EQ(result[3], 4.0);
}

TEST(MatrixTest, Exp)
{
  Matrix<real_t> mat({{0.0, 2.0}, {1.0, 3.0}});
  Matrix<real_t> result = exp(mat);
  EXPECT_DOUBLE_EQ(result[0], std::exp(0.0));
  EXPECT_DOUBLE_EQ(result[1], std::exp(1.0));
  EXPECT_DOUBLE_EQ(result[2], std::exp(2.0));
  EXPECT_DOUBLE_EQ(result[3], std::exp(3.0));
}

TEST(MatrixTest, Log)
{
  Matrix<real_t> mat({{1.0, 3.0}, {2.0, 4.0}});
  Matrix<real_t> result = log(mat);
  EXPECT_DOUBLE_EQ(result[0], std::log(1.0));
  EXPECT_DOUBLE_EQ(result[1], std::log(2.0));
  EXPECT_DOUBLE_EQ(result[2], std::log(3.0));
  EXPECT_DOUBLE_EQ(result[3], std::log(4.0));
}

TEST(MatrixTest, TrigonometricFunctions)
{
  Matrix<real_t> mat({{0.0, M_PI}, {M_PI / 2, 3 * M_PI / 2}});
  Matrix<real_t> sinResult = sin(mat);
  Matrix<real_t> cosResult = cos(mat);

  EXPECT_NEAR(sinResult[0], 0.0, 1e-10);
  EXPECT_NEAR(sinResult[1], 1.0, 1e-10);
  EXPECT_NEAR(cosResult[0], 1.0, 1e-10);
  EXPECT_NEAR(cosResult[1], 0.0, 1e-10);
}

TEST(MatrixTest, Sqrt)
{
  Matrix<real_t> mat({{1.0, 9.0}, {4.0, 16.0}});
  Matrix<real_t> result = sqrt(mat);
  EXPECT_DOUBLE_EQ(result[0], 1.0);
  EXPECT_DOUBLE_EQ(result[1], 2.0);
  EXPECT_DOUBLE_EQ(result[2], 3.0);
  EXPECT_DOUBLE_EQ(result[3], 4.0);
}

TEST(MatrixTest, SqrtNegativeThrowsError)
{
  Matrix<real_t> mat({{-1.0, 9.0}, {4.0, 16.0}});
  // sqrt of negative numbers should throw an error
  EXPECT_THROW(sqrt(mat), MathLibError);
}

TEST(MatrixTest, Elpow)
{
  Matrix<real_t> mat({{2.0, 4.0}, {3.0, 5.0}});
  Matrix<real_t> result = elpow(mat, 2.0);
  EXPECT_DOUBLE_EQ(result[0], 4.0);
  EXPECT_DOUBLE_EQ(result[1], 9.0);
  EXPECT_DOUBLE_EQ(result[2], 16.0);
  EXPECT_DOUBLE_EQ(result[3], 25.0);
}

TEST(MatrixTest, Round)
{
  Matrix<real_t> mat({{1.234, 3.891}, {2.567, 4.123}});
  Matrix<real_t> result = round(mat, 1);
  EXPECT_DOUBLE_EQ(result[0], 1.2);
  EXPECT_DOUBLE_EQ(result[1], 2.6);
  EXPECT_DOUBLE_EQ(result[2], 3.9);
  EXPECT_DOUBLE_EQ(result[3], 4.1);
}

TEST(MatrixTest, Sum)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  real_t result = sum(mat);
  EXPECT_DOUBLE_EQ(result, 10.0);
  Matrix<real_t> rowSum = sum(mat, 1);
  EXPECT_DOUBLE_EQ(rowSum(1), 4.0);
  EXPECT_DOUBLE_EQ(rowSum(2), 6.0);
  Matrix<real_t> colSum = sum(mat, 2);
  EXPECT_DOUBLE_EQ(colSum(1), 3.0);
  EXPECT_DOUBLE_EQ(colSum(2), 7.0);
}

TEST(MatrixTest, Mean)
{
  Matrix<real_t> mat({{1.0, 3.0}, {2.0, 4.0}});
  real_t result = mean(mat);
  EXPECT_DOUBLE_EQ(result, 2.5);
}

TEST(MatrixTest, Max)
{
  Matrix<real_t> mat({{1.0, 3.0}, {5.0, 2.0}});
  EXPECT_DOUBLE_EQ(max(mat), 5.0);
  Matrix<real_t> rowMax = max(mat, 1);
  EXPECT_DOUBLE_EQ(rowMax(1), 5.0);
  EXPECT_DOUBLE_EQ(rowMax(2), 3.0);
  Matrix<real_t> colMax = max(mat, 2);
  EXPECT_DOUBLE_EQ(colMax(1), 3.0);
  EXPECT_DOUBLE_EQ(colMax(2), 5.0);
}

TEST(MatrixTest, Min)
{
  Matrix<real_t> mat({{1.0, 3.0}, {5.0, 2.0}});
  EXPECT_DOUBLE_EQ(min(mat), 1.0);
  Matrix<real_t> rowMin = min(mat, 1);
  EXPECT_DOUBLE_EQ(rowMin(1), 1.0);
  EXPECT_DOUBLE_EQ(rowMin(2), 2.0);
  Matrix<real_t> colMin = min(mat, 2);
  EXPECT_DOUBLE_EQ(colMin(1), 1.0);
  EXPECT_DOUBLE_EQ(colMin(2), 2.0);
}

TEST(MatrixTest, Elmul)
{
  Matrix<real_t> mat1({{1.0, 3.0}, {2.0, 4.0}});
  Matrix<real_t> mat2({{2.0, 4.0}, {3.0, 5.0}});
  Matrix<real_t> result = elmul(mat1, mat2);
  EXPECT_DOUBLE_EQ(result[0], 2.0);
  EXPECT_DOUBLE_EQ(result[1], 6.0);
  EXPECT_DOUBLE_EQ(result[2], 12.0);
  EXPECT_DOUBLE_EQ(result[3], 20.0);
}

TEST(MatrixTest, Eldiv)
{
  Matrix<real_t> mat1({{10.0, 30.0}, {20.0, 40.0}});
  Matrix<real_t> mat2({{2.0, 5.0}, {4.0, 8.0}});
  Matrix<real_t> result = eldiv(mat1, mat2);
  EXPECT_DOUBLE_EQ(result[0], 5.0);
  EXPECT_DOUBLE_EQ(result[1], 5.0);
  EXPECT_DOUBLE_EQ(result[2], 6.0);
  EXPECT_DOUBLE_EQ(result[3], 5.0);
}

TEST(MatrixTest, AllAny)
{
  Matrix<bool> mat1({{true, true}, {true, true}});
  Matrix<bool> mat2({{true, false}, {true, false}});
  Matrix<bool> mat3({{false, false}, {false, false}});

  EXPECT_TRUE(all(mat1));
  EXPECT_FALSE(all(mat2));
  EXPECT_FALSE(all(mat3));

  EXPECT_TRUE(any(mat1));
  EXPECT_TRUE(any(mat2));
  EXPECT_FALSE(any(mat3));
}

// TODO fix test and also test for exceptions of invalid sizes

TEST(MatrixTest, Hcat)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0}, {7.0, 8.0}});
  Matrix<real_t> result = hcat(mat1, mat2);
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 4);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(1, 3), 5.0);
  EXPECT_DOUBLE_EQ(result(1, 4), 6.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(2, 3), 7.0);
  EXPECT_DOUBLE_EQ(result(2, 4), 8.0);
  Matrix<real_t> mat3({{1.0, 2.0, 3.0}});
  EXPECT_THROW(hcat(mat1, mat3), IncompatibleSizeError);
}

TEST(MatrixTest, Vcat)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0}, {7.0, 8.0}});
  Matrix<real_t> result = vcat(mat1, mat2);
  EXPECT_EQ(result.rows(), 4);
  EXPECT_EQ(result.cols(), 2);
  EXPECT_DOUBLE_EQ(result(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(result(1, 2), 2.0);
  EXPECT_DOUBLE_EQ(result(2, 1), 3.0);
  EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
  EXPECT_DOUBLE_EQ(result(3, 1), 5.0);
  EXPECT_DOUBLE_EQ(result(3, 2), 6.0);
  EXPECT_DOUBLE_EQ(result(4, 1), 7.0);
  EXPECT_DOUBLE_EQ(result(4, 2), 8.0);
  Matrix<real_t> mat3({{1.0, 2.0, 3.0}});
  EXPECT_THROW(vcat(mat1, mat3), IncompatibleSizeError);
}

TEST(MatrixTest, Cat)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0}, {7.0, 8.0}});
  Matrix<real_t> resultV = cat(mat1, mat2, 1);
  EXPECT_TRUE(all(resultV == vcat(mat1, mat2)));
  Matrix<real_t> resultH = cat(mat1, mat2, 2);
  EXPECT_TRUE(all(resultH == hcat(mat1, mat2)));
}

TEST(MatrixTest, FuzzyComparisons)
{
  Matrix<real_t> mat1({{1, 2}, {3, 4}});
  Matrix<real_t> mat2({{1, 3}, {2, -4}});

  Matrix<bool> eq = isFuzzyEqual(mat1, mat2);
  EXPECT_TRUE(eq(1, 1));
  EXPECT_FALSE(eq(1, 2));
  EXPECT_FALSE(eq(2, 1));
  EXPECT_FALSE(eq(2, 2));

  Matrix<bool> gt = isStrictFuzzyGreater(mat1, mat2);
  EXPECT_FALSE(gt(1, 1));
  EXPECT_FALSE(gt(1, 2));
  EXPECT_TRUE(gt(2, 1));
  EXPECT_TRUE(gt(2, 2));

  Matrix<bool> gteq = isFuzzyGreater(mat1, mat2);
  EXPECT_FALSE(gteq(1, 1));
  EXPECT_FALSE(gteq(1, 2));
  EXPECT_TRUE(gteq(2, 1));
  EXPECT_TRUE(gteq(2, 2));

  Matrix<bool> lt = isStrictFuzzySmaller(mat1, mat2);
  EXPECT_FALSE(lt(1, 1));
  EXPECT_TRUE(lt(1, 2));
  EXPECT_FALSE(lt(2, 1));
  EXPECT_FALSE(lt(2, 2));

  Matrix<bool> lteq = isFuzzySmaller(mat1, mat2);
  EXPECT_FALSE(lteq(1, 1));
  EXPECT_TRUE(lteq(1, 2));
  EXPECT_FALSE(lteq(2, 1));
  EXPECT_FALSE(lteq(2, 2));
}

TEST(MatrixTest, IncompatibleSizes)
{
  Matrix<real_t> mat1(2, 2);
  Matrix<real_t> mat2(3, 3);
  EXPECT_THROW(mat1 + mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 - mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 * mat2, IncompatibleSizeError);
}

TEST(MatrixTest, IndexOutOfRange)
{
  Matrix<real_t> mat(2, 2);
  EXPECT_THROW(mat[4], IndexOutOfRangeError);
  EXPECT_THROW(mat(0, 0), IndexOutOfRangeError);
  EXPECT_THROW(mat(0, 1), IndexOutOfRangeError);
  EXPECT_THROW(mat(1, 0), IndexOutOfRangeError);
  EXPECT_THROW(mat(3, 3), IndexOutOfRangeError);
}

TEST(MatrixTest, MixedTypeOperations)
{
  Matrix<int64_t> mat1({{1, 2}, {3, 4}});
  Matrix<real_t> mat2({{1.5, 2.5}, {3.5, 4.5}});
  auto result = mat1 + mat2;
  EXPECT_DOUBLE_EQ(result(1, 1), 2.5);
  EXPECT_DOUBLE_EQ(result(1, 2), 4.5);
  EXPECT_DOUBLE_EQ(result(2, 1), 6.5);
  EXPECT_DOUBLE_EQ(result(2, 2), 8.5);
}

TEST(MatrixTest, Norm)
{
  Matrix<real_t> I = eye(2);
  EXPECT_DOUBLE_EQ(norm(I), std::sqrt(2));
  Matrix<real_t> A{{1.0, 2.0}, {3.0, 4.0}};
  EXPECT_DOUBLE_EQ(norm(A), std::sqrt(trace(transpose(A) * A)));
}

TEST(MatrixTest, Vecnorm)
{
  // Test matrix for column norms (dim=1)
  Matrix<real_t> A({{3.0, 0.0}, {4.0, 5.0}, {0.0, 12.0}});

  // Test 2-norm along columns (default)
  Matrix<real_t> colNorms2 = vecnorm(A, 2.0, 1);
  EXPECT_EQ(colNorms2.rows(), 1);
  EXPECT_EQ(colNorms2.cols(), 2);
  EXPECT_DOUBLE_EQ(colNorms2(1, 1), 5.0);   // sqrt(3^2 + 4^2 + 0^2) = sqrt(25) = 5
  EXPECT_DOUBLE_EQ(colNorms2(1, 2), 13.0);  // sqrt(0^2 + 5^2 + 12^2) = sqrt(169) = 13

  // Test 1-norm along columns
  Matrix<real_t> colNorms1 = vecnorm(A, 1.0, 1);
  EXPECT_EQ(colNorms1.rows(), 1);
  EXPECT_EQ(colNorms1.cols(), 2);
  EXPECT_DOUBLE_EQ(colNorms1(1, 1), 7.0);   // |3| + |4| + |0| = 7
  EXPECT_DOUBLE_EQ(colNorms1(1, 2), 17.0);  // |0| + |5| + |12| = 17

  // Test infinity norm along columns (max absolute value)
  Matrix<real_t> colNormsInf = vecnorm(A, std::numeric_limits<real_t>::infinity(), 1);
  EXPECT_EQ(colNormsInf.rows(), 1);
  EXPECT_EQ(colNormsInf.cols(), 2);
  EXPECT_DOUBLE_EQ(colNormsInf(1, 1), 4.0);   // max(|3|, |4|, |0|) = 4
  EXPECT_DOUBLE_EQ(colNormsInf(1, 2), 12.0);  // max(|0|, |5|, |12|) = 12

  // Test -infinity norm along columns (min absolute value)
  Matrix<real_t> colNormsNegInf = vecnorm(A, -std::numeric_limits<real_t>::infinity(), 1);
  EXPECT_EQ(colNormsNegInf.rows(), 1);
  EXPECT_EQ(colNormsNegInf.cols(), 2);
  EXPECT_DOUBLE_EQ(colNormsNegInf(1, 1), 0.0);  // min(|3|, |4|, |0|) = 0
  EXPECT_DOUBLE_EQ(colNormsNegInf(1, 2), 0.0);  // min(|0|, |5|, |12|) = 0

  // Test 2-norm along rows (dim=2)
  Matrix<real_t> B({{1.0, 2.0, 2.0}, {3.0, 4.0, 0.0}});
  Matrix<real_t> rowNorms2 = vecnorm(B, 2.0, 2);
  EXPECT_EQ(rowNorms2.rows(), 2);
  EXPECT_EQ(rowNorms2.cols(), 1);
  EXPECT_DOUBLE_EQ(rowNorms2(1, 1), 3.0);  // sqrt(1^2 + 2^2 + 2^2) = sqrt(9) = 3
  EXPECT_DOUBLE_EQ(rowNorms2(2, 1), 5.0);  // sqrt(3^2 + 4^2 + 0^2) = sqrt(25) = 5

  // Test 1-norm along rows
  Matrix<real_t> rowNorms1 = vecnorm(B, 1.0, 2);
  EXPECT_EQ(rowNorms1.rows(), 2);
  EXPECT_EQ(rowNorms1.cols(), 1);
  EXPECT_DOUBLE_EQ(rowNorms1(1, 1), 5.0);  // |1| + |2| + |2| = 5
  EXPECT_DOUBLE_EQ(rowNorms1(2, 1), 7.0);  // |3| + |4| + |0| = 7

  // Test infinity norm along rows
  Matrix<real_t> rowNormsInf = vecnorm(B, std::numeric_limits<real_t>::infinity(), 2);
  EXPECT_EQ(rowNormsInf.rows(), 2);
  EXPECT_EQ(rowNormsInf.cols(), 1);
  EXPECT_DOUBLE_EQ(rowNormsInf(1, 1), 2.0);  // max(|1|, |2|, |2|) = 2
  EXPECT_DOUBLE_EQ(rowNormsInf(2, 1), 4.0);  // max(|3|, |4|, |0|) = 4

  // Test p-norm with p=3 along columns
  Matrix<real_t> C({{1.0, 2.0}, {2.0, 1.0}, {3.0, 0.0}});
  Matrix<real_t> colNorms3 = vecnorm(C, 3.0, 1);
  EXPECT_EQ(colNorms3.rows(), 1);
  EXPECT_EQ(colNorms3.cols(), 2);
  // Column 1: (1^3 + 2^3 + 3^3)^(1/3) = (1 + 8 + 27)^(1/3) = 36^(1/3)
  EXPECT_NEAR(colNorms3(1, 1), std::pow(36.0, 1.0 / 3.0), 1e-10);
  // Column 2: (2^3 + 1^3 + 0^3)^(1/3) = (8 + 1 + 0)^(1/3) = 9^(1/3)
  EXPECT_NEAR(colNorms3(1, 2), std::pow(9.0, 1.0 / 3.0), 1e-10);

  // Test with negative values
  Matrix<real_t> D({{-3.0, 4.0}, {-5.0, -12.0}});
  Matrix<real_t> DNorms = vecnorm(D, 2.0, 1);
  EXPECT_NEAR(DNorms(1, 1), std::sqrt(9.0 + 25.0), 1e-10);    // sqrt(34)
  EXPECT_NEAR(DNorms(1, 2), std::sqrt(16.0 + 144.0), 1e-10);  // sqrt(160) = 4*sqrt(10)

  // Test zero matrix
  Matrix<real_t> Z(2, 3);
  Matrix<real_t> ZNorms1 = vecnorm(Z, 2.0, 1);
  Matrix<real_t> ZNorms2 = vecnorm(Z, 2.0, 2);
  EXPECT_EQ(ZNorms1.rows(), 1);
  EXPECT_EQ(ZNorms1.cols(), 3);
  EXPECT_EQ(ZNorms2.rows(), 2);
  EXPECT_EQ(ZNorms2.cols(), 1);
  for (size_t i = 0; i < ZNorms1.length(); ++i) {
    EXPECT_DOUBLE_EQ(ZNorms1[i], 0.0);
  }
  for (size_t i = 0; i < ZNorms2.length(); ++i) {
    EXPECT_DOUBLE_EQ(ZNorms2[i], 0.0);
  }

  // Test that p <= 0 throws error
  Matrix<real_t> E(2, 2);
  EXPECT_THROW(vecnorm(E, 0.0), MathLibError);
  EXPECT_THROW(vecnorm(E, -1.0), MathLibError);
}

TEST(MatrixTest, QRDecomposition_Square)
{
  // Test with a square matrix
  Matrix<real_t> A({{12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}});

  auto qr_result = qr(A);
  Matrix<real_t> Q = qr_result.Q;
  Matrix<real_t> R = qr_result.R;

  // Verify dimensions
  EXPECT_EQ(Q.rows(), 3);
  EXPECT_EQ(Q.cols(), 3);
  EXPECT_EQ(R.rows(), 3);
  EXPECT_EQ(R.cols(), 3);

  // Verify Q*R = A
  Matrix<real_t> QR = Q * R;
  for (size_t i = 0; i < A.length(); ++i) {
    EXPECT_NEAR(QR[i], A[i], 1e-10);
  }

  // Verify Q is orthogonal (Q'*Q = I)
  Matrix<real_t> QtQ = transpose(Q) * Q;
  Matrix<real_t> I = eye<real_t>(3);
  for (size_t i = 0; i < I.length(); ++i) {
    EXPECT_NEAR(QtQ[i], I[i], 1e-10);
  }

  // Verify R is upper triangular
  EXPECT_NEAR(R(2, 1), 0.0, 1e-10);
  EXPECT_NEAR(R(3, 1), 0.0, 1e-10);
  EXPECT_NEAR(R(3, 2), 0.0, 1e-10);
}

TEST(MatrixTest, QRDecomposition_Rectangular)
{
  // Test with a tall rectangular matrix (more rows than columns)
  // Use a matrix with linearly independent columns
  Matrix<real_t> A({{1.0, 0.0, 2.0}, {0.0, 3.0, 1.0}, {0.0, 4.0, 2.0}, {1.0, 1.0, 0.0}});

  auto qr_result = qr(A);
  Matrix<real_t> Q = qr_result.Q;
  Matrix<real_t> R = qr_result.R;

  // Verify dimensions
  EXPECT_EQ(Q.rows(), 4);
  EXPECT_EQ(Q.cols(), 3);
  EXPECT_EQ(R.rows(), 3);
  EXPECT_EQ(R.cols(), 3);

  // Verify Q*R = A
  Matrix<real_t> QR = Q * R;
  for (size_t i = 0; i < A.length(); ++i) {
    EXPECT_NEAR(QR[i], A[i], 1e-10);
  }

  // Verify Q has orthonormal columns (Q'*Q = I)
  Matrix<real_t> QtQ = transpose(Q) * Q;
  Matrix<real_t> I = eye<real_t>(3);
  for (size_t i = 0; i < I.length(); ++i) {
    EXPECT_NEAR(QtQ[i], I[i], 1e-10);
  }
}

TEST(MatrixTest, QRDecomposition_Identity)
{
  // QR of identity should be Q=I, R=I
  Matrix<real_t> I = eye<real_t>(3);

  auto qr_result = qr(I);
  Matrix<real_t> Q = qr_result.Q;
  Matrix<real_t> R = qr_result.R;

  // Q should be identity
  for (size_t i = 0; i < I.length(); ++i) {
    EXPECT_NEAR(Q[i], I[i], 1e-10);
  }

  // R should be identity
  for (size_t i = 0; i < I.length(); ++i) {
    EXPECT_NEAR(R[i], I[i], 1e-10);
  }
}

TEST(MatrixTest, QRDecomposition_SimpleMatrix)
{
  // Simple 2x2 matrix
  Matrix<real_t> A({{3.0, 0.0}, {4.0, 5.0}});

  auto qr_result = qr(A);
  Matrix<real_t> Q = qr_result.Q;
  Matrix<real_t> R = qr_result.R;

  // Verify Q*R = A
  Matrix<real_t> QR = Q * R;
  for (size_t i = 0; i < A.length(); ++i) {
    EXPECT_NEAR(QR[i], A[i], 1e-10);
  }

  // Verify R is upper triangular
  EXPECT_NEAR(R(2, 1), 0.0, 1e-10);
}

TEST(MatrixTest, QRDecomposition_RankDeficient)
{
  // Rank deficient matrix (column 2 is zero)
  Matrix<real_t> A({{1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, {4.0, 5.0, 6.0}});

  // Should throw error for rank deficient matrix by default
  EXPECT_THROW(qr(A), MathLibError);

  // Should succeed with checkRank=false
  EXPECT_NO_THROW(qr(A, false));
}

TEST(MatrixTest, Orth)
{
  // Create a matrix with linearly dependent columns
  Matrix<real_t> A({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}});  // Column 2 = Column 0 + Column 1

  Matrix<real_t> Q = orth(A);

  // Should only have 2 orthonormal columns (rank = 2)
  EXPECT_EQ(Q.rows(), 3);
  EXPECT_EQ(Q.cols(), 2);

  // Verify columns are orthonormal
  Matrix<real_t> QtQ = transpose(Q) * Q;
  Matrix<real_t> I = eye<real_t>(2);

  for (size_t i = 0; i < I.length(); ++i) {
    EXPECT_NEAR(QtQ[i], I[i], 1e-10);
  }
}

TEST(MatrixTest, Rank_FullRank)
{
  // Full rank 3x3 matrix
  Matrix<real_t> A({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});

  EXPECT_EQ(rank(A), 3);
  EXPECT_TRUE(isFullRank(A));
}

TEST(MatrixTest, Rank_RankDeficient)
{
  // Rank deficient matrix (column 2 is zero)
  Matrix<real_t> A({{1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, {4.0, 5.0, 6.0}});

  EXPECT_EQ(rank(A), 2);
  EXPECT_FALSE(isFullRank(A));
}

TEST(MatrixTest, Rank_Rectangular)
{
  // 4x3 matrix with rank 3
  Matrix<real_t> A({{1.0, 0.0, 2.0}, {0.0, 3.0, 1.0}, {0.0, 4.0, 2.0}, {1.0, 1.0, 0.0}});

  EXPECT_EQ(rank(A), 3);
  EXPECT_TRUE(isFullRank(A));
}

TEST(MatrixTest, Rank_ZeroMatrix)
{
  Matrix<real_t> A(3, 3);  // Zero matrix

  EXPECT_EQ(rank(A), 0);
  EXPECT_FALSE(isFullRank(A));
}

TEST(MatrixTest, Dot)
{
  Matrix<real_t> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<real_t> B = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

  Matrix<real_t> result1 = dot(A, B, 1);
  EXPECT_EQ(result1.rows(), 1);
  EXPECT_EQ(result1.cols(), 3);
  EXPECT_DOUBLE_EQ(result1(1), 54.0);
  EXPECT_DOUBLE_EQ(result1(2), 57.0);
  EXPECT_DOUBLE_EQ(result1(3), 54.0);

  Matrix<real_t> result2 = dot(A, B, 2);
  EXPECT_EQ(result2.rows(), 3);
  EXPECT_EQ(result2.cols(), 1);
  EXPECT_DOUBLE_EQ(result2(1), 46.0);
  EXPECT_DOUBLE_EQ(result2(2), 73.0);
  EXPECT_DOUBLE_EQ(result2(3), 46.0);

  Matrix<real_t> C = eye(2);
  EXPECT_THROW(dot(A, C), IncompatibleSizeError);
  EXPECT_THROW(dot(C, A), IncompatibleSizeError);
}

TEST(MatrixTest, ColDot)
{
  Matrix<real_t> A({{1.0, 2.0}, {3.0, 4.0}});
  EXPECT_DOUBLE_EQ(colDot(A, 1, 1), 10.0);
  EXPECT_DOUBLE_EQ(colDot(A, 1, 2), 14.0);
  EXPECT_DOUBLE_EQ(colDot(A, 2, 1), 14.0);
  EXPECT_DOUBLE_EQ(colDot(A, 2, 2), 20.0);
  EXPECT_THROW(colDot(A, 3, 1), IndexOutOfRangeError);
  EXPECT_THROW(colDot(A, 1, 3), IndexOutOfRangeError);
  EXPECT_THROW(colDot(A, 3, 3), IndexOutOfRangeError);
}

TEST(MatrixTest, RowDot)
{
  Matrix<real_t> A({{1.0, 2.0}, {3.0, 4.0}});
  EXPECT_DOUBLE_EQ(rowDot(A, 1, 1), 5.0);
  EXPECT_DOUBLE_EQ(rowDot(A, 1, 2), 11.0);
  EXPECT_DOUBLE_EQ(rowDot(A, 2, 1), 11.0);
  EXPECT_DOUBLE_EQ(rowDot(A, 2, 2), 25.0);
  EXPECT_THROW(rowDot(A, 3, 1), IndexOutOfRangeError);
  EXPECT_THROW(rowDot(A, 1, 3), IndexOutOfRangeError);
  EXPECT_THROW(rowDot(A, 3, 3), IndexOutOfRangeError);
}

TEST(MatrixTest, Normalize)
{
  Matrix<real_t> A({{3.0, 4.0}, {0.0, 5.0}, {0.0, 12.0}});
  Matrix<real_t> NA = normalize(A, 2, 1);

  // Each column should have norm 1
  Matrix<real_t> normsA = vecnorm(NA, 2, 1);
  for (size_t i = 0; i < normsA.length(); ++i) {
    EXPECT_NEAR(normsA[i], 1.0, 1e-10);
  }

  Matrix<real_t> B({{3.0, 4.0, 0.0}, {5.0, 0.0, 12.0}});
  Matrix<real_t> NB = normalize(B, 2, 2);

  // Each row should have norm 1
  Matrix<real_t> normsB = vecnorm(NB, 2, 2);
  for (size_t i = 0; i < normsB.length(); ++i) {
    EXPECT_NEAR(normsB[i], 1.0, 1e-10);
  }

  Matrix<real_t> Z(2, 2);  // Zero matrix
  EXPECT_THROW(normalize(Z, 2, 1), MathLibError);
}

TEST(MatrixTest, ConstructorEdgeCases)
{
  // Empty initializer list should create an empty matrix
  std::initializer_list<std::initializer_list<real_t>> empty_list = {};
  Matrix<real_t> mat(empty_list);
  EXPECT_EQ(mat.rows(), 0);
  EXPECT_EQ(mat.cols(), 0);
  EXPECT_TRUE(mat.isEmpty());
}

TEST(MatrixTest, InitializerListMismatchedRowSizes)
{
  // Rows with different column counts should throw an error
  EXPECT_THROW((Matrix<real_t>({{1.0, 2.0}, {3.0, 4.0, 5.0}})), MathLibError);
}

TEST(MatrixTest, EmptyVectorListConstructor)
{
  // Empty vector list should create an empty matrix
  std::initializer_list<Vector<real_t>> empty_list;
  Matrix<real_t> mat(empty_list);
  EXPECT_EQ(mat.rows(), 0);
  EXPECT_EQ(mat.cols(), 0);
  EXPECT_TRUE(mat.isEmpty());
}

TEST(MatrixTest, VectorListMismatchedLengths)
{
  // Vectors with different lengths should throw an error
  Vector<real_t> v1 = {1.0, 2.0, 3.0};
  Vector<real_t> v2 = {4.0, 5.0};
  EXPECT_THROW((Matrix<real_t>({v1, v2})), MathLibError);
}

TEST(MatrixTest, EmptyVectorConstructor)
{
  // Constructing from an empty vector
  Vector<real_t> empty_vec;
  Matrix<real_t> mat(empty_vec);
  EXPECT_EQ(mat.rows(), 0);
  EXPECT_EQ(mat.cols(), 0);
  EXPECT_TRUE(mat.isEmpty());
}

TEST(MatrixTest, SizeConstructorWithValue)
{
  // Constructor with size and fill value
  Matrix<real_t> mat(2, 3, 7.5);
  EXPECT_EQ(mat.rows(), 2);
  EXPECT_EQ(mat.cols(), 3);
  for (size_t i = 0; i < mat.length(); ++i) {
    EXPECT_DOUBLE_EQ(mat[i], 7.5);
  }
}

TEST(MatrixTest, DivisionByZero)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});
  EXPECT_THROW(mat / 0.0, DivByZeroError);
  EXPECT_THROW(mat / static_cast<int64_t>(0), DivByZeroError);
}

TEST(MatrixTest, EldivDivisionByZero)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{1.0, 0.0}, {2.0, 3.0}});
  EXPECT_THROW(eldiv(mat1, mat2), DivByZeroError);
}

TEST(MatrixTest, NegativeIndexingEdgeCases)
{
  Matrix<real_t> mat({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});

  // Test negative indexing at boundaries
  EXPECT_DOUBLE_EQ(mat(-1, -1), 6.0);  // Last element
  EXPECT_DOUBLE_EQ(mat(-2, -3), 1.0);  // First element

  // Test out of range negative indices
  EXPECT_THROW(mat(-3, 1), IndexOutOfRangeError);
  EXPECT_THROW(mat(1, -4), IndexOutOfRangeError);
}

TEST(MatrixTest, UnaryOperatorsOnEmpty)
{
  Matrix<real_t> empty;
  Matrix<real_t> plus_result = +empty;
  Matrix<real_t> minus_result = -empty;

  EXPECT_TRUE(plus_result.isEmpty());
  EXPECT_TRUE(minus_result.isEmpty());
}

TEST(MatrixTest, ComparisonIncompatibleSizes)
{
  Matrix<real_t> mat1(2, 2);
  Matrix<real_t> mat2(2, 3);

  EXPECT_THROW(mat1 == mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 != mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 < mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 > mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 <= mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 >= mat2, IncompatibleSizeError);
}

TEST(MatrixTest, LogicalOperatorsIncompatibleSizes)
{
  Matrix<bool> mat1(2, 2);
  Matrix<bool> mat2(2, 3);

  EXPECT_THROW(mat1 && mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 || mat2, IncompatibleSizeError);
  EXPECT_THROW(mat1 ^ mat2, IncompatibleSizeError);
}

TEST(MatrixTest, MixedTypeAdditionSubtraction)
{
  Matrix<int64_t> mat_int({{1, 2}, {3, 4}});
  Matrix<real_t> mat_real({{1.5, 2.5}, {3.5, 4.5}});

  // int + real should return real
  auto result1 = mat_int + mat_real;
  EXPECT_DOUBLE_EQ(result1(1, 1), 2.5);
  EXPECT_DOUBLE_EQ(result1(2, 2), 8.5);

  // real - int should return real
  auto result2 = mat_real - mat_int;
  EXPECT_DOUBLE_EQ(result2(1, 1), 0.5);
  EXPECT_DOUBLE_EQ(result2(2, 2), 0.5);
}

TEST(MatrixTest, MixedTypeMultiplication)
{
  Matrix<int64_t> mat_int({{1, 2}, {3, 4}});
  Matrix<real_t> mat_real({{1.5, 2.5}, {3.5, 4.5}});

  // int * real matrix multiplication
  auto result = mat_int * mat_real;
  EXPECT_EQ(result.rows(), 2);
  EXPECT_EQ(result.cols(), 2);

  // Verify result values: [1 2; 3 4] * [1.5 2.5; 3.5 4.5]
  // First row, first col: 1*1.5 + 2*3.5 = 1.5 + 7.0 = 8.5
  EXPECT_DOUBLE_EQ(result(1, 1), 8.5);
}

TEST(MatrixTest, MixedTypeVectorMultiplication)
{
  Matrix<int64_t> mat_int({{1, 2}, {3, 4}});
  Vector<real_t> vec_real = {1.5, 2.5};

  auto result = mat_int * vec_real;
  EXPECT_EQ(result.length(), 2);
  // [1 2; 3 4] * [1.5; 2.5] = [6.5; 14.5]
  EXPECT_DOUBLE_EQ(result[0], 6.5);
  EXPECT_DOUBLE_EQ(result[1], 14.5);
}

TEST(MatrixTest, VectorMatrixMultiplicationEdgeCases)
{
  // Vector * Matrix treats vector as column vector
  // So a 3-element vector times a 1x2 matrix = 3x1 * 1x2 = 3x2 matrix
  Vector<real_t> vec = {1.0, 2.0, 3.0};
  Matrix<real_t> mat({{4.0, 5.0}});  // 1x2 matrix
  auto result = vec * mat;
  EXPECT_EQ(result.rows(), 3);
  EXPECT_EQ(result.cols(), 2);

  // Incompatible sizes: vector length doesn't match matrix rows
  Vector<real_t> vec_bad = {1.0, 2.0};
  Matrix<real_t> mat_bad({{1.0, 2.0}, {3.0, 4.0}});  // 2x2 matrix
  EXPECT_THROW(vec_bad * mat_bad, IncompatibleSizeError);
}

TEST(MatrixTest, CatInvalidDimension)
{
  Matrix<real_t> mat1({{1.0, 2.0}, {3.0, 4.0}});
  Matrix<real_t> mat2({{5.0, 6.0}, {7.0, 8.0}});

  // Invalid dimension (not 1 or 2)
  EXPECT_THROW(cat(mat1, mat2, 0), MathLibError);
  EXPECT_THROW(cat(mat1, mat2, 3), MathLibError);
}

TEST(MatrixTest, HcatIncompatibleRows)
{
  Matrix<real_t> mat1(2, 3);
  Matrix<real_t> mat2(3, 2);

  EXPECT_THROW(hcat(mat1, mat2), IncompatibleSizeError);
}

TEST(MatrixTest, VcatIncompatibleCols)
{
  Matrix<real_t> mat1(2, 3);
  Matrix<real_t> mat2(2, 2);

  EXPECT_THROW(vcat(mat1, mat2), IncompatibleSizeError);
}

TEST(MatrixTest, MultiplicationIncompatibleDimensions)
{
  Matrix<real_t> mat1(2, 3);
  Matrix<real_t> mat2(2, 2);

  // 2x3 * 2x2 should fail (inner dimensions don't match)
  EXPECT_THROW(mat1 * mat2, IncompatibleSizeError);
}

TEST(MatrixTest, MatrixVectorMultiplicationIncompatible)
{
  Matrix<real_t> mat(2, 3);
  Vector<real_t> vec(2);  // Vector length doesn't match matrix cols

  EXPECT_THROW(mat * vec, IncompatibleSizeError);
}

TEST(MatrixTest, SumInvalidDimension)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});

  EXPECT_THROW(sum(mat, 0), MathLibError);
  EXPECT_THROW(sum(mat, 3), MathLibError);
}

TEST(MatrixTest, MaxMinInvalidDimension)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});

  EXPECT_THROW(max(mat, 0), MathLibError);
  EXPECT_THROW(max(mat, 3), MathLibError);
  EXPECT_THROW(min(mat, 0), MathLibError);
  EXPECT_THROW(min(mat, 3), MathLibError);
}

TEST(MatrixTest, ElmulIncompatibleSizes)
{
  Matrix<real_t> mat1(2, 2);
  Matrix<real_t> mat2(2, 3);

  EXPECT_THROW(elmul(mat1, mat2), IncompatibleSizeError);
}

TEST(MatrixTest, EldivIncompatibleSizes)
{
  Matrix<real_t> mat1(2, 2);
  Matrix<real_t> mat2(2, 3);

  EXPECT_THROW(eldiv(mat1, mat2), IncompatibleSizeError);
}

TEST(MatrixTest, ScalarMatrixOperations)
{
  Matrix<real_t> mat({{1.0, 2.0}, {3.0, 4.0}});

  // Scalar * Matrix
  auto result1 = 2.0 * mat;
  EXPECT_DOUBLE_EQ(result1(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(result1(2, 2), 8.0);

  // Matrix * Scalar
  auto result2 = mat * 3.0;
  EXPECT_DOUBLE_EQ(result2(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(result2(2, 2), 12.0);

  // Scalar + Matrix
  auto result3 = 5.0 + mat;
  EXPECT_DOUBLE_EQ(result3(1, 1), 6.0);
  EXPECT_DOUBLE_EQ(result3(2, 2), 9.0);

  // Matrix + Scalar
  auto result4 = mat + 10.0;
  EXPECT_DOUBLE_EQ(result4(1, 1), 11.0);
  EXPECT_DOUBLE_EQ(result4(2, 2), 14.0);

  // Scalar - Matrix
  auto result5 = 10.0 - mat;
  EXPECT_DOUBLE_EQ(result5(1, 1), 9.0);
  EXPECT_DOUBLE_EQ(result5(2, 2), 6.0);

  // Matrix - Scalar
  auto result6 = mat - 1.0;
  EXPECT_DOUBLE_EQ(result6(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(result6(2, 2), 3.0);
}

TEST(MatrixTest, IntScalarRealMatrixOperations)
{
  Matrix<real_t> mat({{1.5, 2.5}, {3.5, 4.5}});

  // int * real Matrix
  auto result1 = 2 * mat;
  EXPECT_DOUBLE_EQ(result1(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(result1(2, 2), 9.0);

  // real Matrix * int
  auto result2 = mat * 2;
  EXPECT_DOUBLE_EQ(result2(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(result2(2, 2), 9.0);
}