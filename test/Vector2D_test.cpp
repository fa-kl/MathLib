#include "DivByZeroError.hpp"
#include "IndexOutOfRangeError.hpp"
#include "MathLibError.hpp"
#include "Vector2D.hpp"

#include <complex>
#include <gtest/gtest.h>

using namespace mathlib;

TEST(Vector2DTest, DefaultConstructor) {
  const Vector2D vec;
  EXPECT_DOUBLE_EQ(vec.x(), 0.0);
  EXPECT_DOUBLE_EQ(vec.y(), 0.0);
  EXPECT_DOUBLE_EQ(x(vec), 0.0);
  EXPECT_DOUBLE_EQ(y(vec), 0.0);
}

TEST(Vector2DTest, Constructor) {
  const Vector2D vec1(1.0, 2.0);
  EXPECT_DOUBLE_EQ(vec1.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec1.y(), 2.0);
  const Vector2D vec2(-2.0, -1.0);
  EXPECT_DOUBLE_EQ(vec2.x(), -2.0);
  EXPECT_DOUBLE_EQ(vec2.y(), -1.0);
}

TEST(Vector2DTest, CopyConstructor) {
  const Vector2D vec1(1.0, 2.0);
  const Vector2D vec2(vec1);
  EXPECT_DOUBLE_EQ(vec2.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 2.0);
}

TEST(Vector2DTest, Copy) {
  const Vector2D vec1(1.0, 2.0);
  Vector2D vec2;
  vec2 = vec1;
  EXPECT_DOUBLE_EQ(vec2.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 2.0);
}

TEST(Vector2DTest, GetterMehtods) {
  Vector2D vec1(1.0, 2.0);
  EXPECT_DOUBLE_EQ(vec1.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec1.y(), 2.0);

  vec1.x() = 0;
  vec1.y() = 0;
  EXPECT_DOUBLE_EQ(vec1[0], 0.0);
  EXPECT_DOUBLE_EQ(vec1[1], 0.0);

  const Vector2D vec2(3.0, 4.0);
  EXPECT_DOUBLE_EQ(vec2.x(), 3.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 4.0);
}

TEST(Vector2DTest, GetterSquareBracketOperator) {
  Vector2D vec1(1.0, 2.0);
  EXPECT_DOUBLE_EQ(vec1[0], 1.0);
  EXPECT_DOUBLE_EQ(vec1[1], 2.0);
  EXPECT_THROW(vec1[2], mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec1(0), mathlib::IndexOutOfRangeError);

  vec1[0] = 0;
  vec1[1] = 0;
  EXPECT_DOUBLE_EQ(vec1[0], 0.0);
  EXPECT_DOUBLE_EQ(vec1[1], 0.0);

  const Vector2D vec2(3.0, 4.0);
  EXPECT_DOUBLE_EQ(vec2[0], 3.0);
  EXPECT_DOUBLE_EQ(vec2[1], 4.0);
  EXPECT_THROW(vec2[2], mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec2(0), mathlib::IndexOutOfRangeError);
}

TEST(Vector2DTest, GetterRoundBracketOperator) {
  Vector2D vec1(1.0, 2.0);
  EXPECT_DOUBLE_EQ(vec1(1), 1.0);
  EXPECT_DOUBLE_EQ(vec1(2), 2.0);
  EXPECT_DOUBLE_EQ(vec1(-2), 1.0);
  EXPECT_DOUBLE_EQ(vec1(-1), 2.0);
  EXPECT_THROW(vec1(0), mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec1(3), mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec1(-3), mathlib::IndexOutOfRangeError);

  vec1(1) = 0;
  vec1(2) = 0;
  EXPECT_DOUBLE_EQ(vec1(1), 0.0);
  EXPECT_DOUBLE_EQ(vec1(2), 0.0);
  vec1(-2) = 2;
  vec1(-1) = 1;
  EXPECT_DOUBLE_EQ(vec1(1), 2.0);
  EXPECT_DOUBLE_EQ(vec1(2), 1.0);

  const Vector2D vec2(3.0, 4.0);
  EXPECT_DOUBLE_EQ(vec2(1), 3.0);
  EXPECT_DOUBLE_EQ(vec2(2), 4.0);
  EXPECT_DOUBLE_EQ(vec2(-2), 3.0);
  EXPECT_DOUBLE_EQ(vec2(-1), 4.0);
  EXPECT_THROW(vec2(0), mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec2(3), mathlib::IndexOutOfRangeError);
  EXPECT_THROW(vec2(-3), mathlib::IndexOutOfRangeError);
}

TEST(Vector2DTest, Norm) {
  const Vector2D zero(0.0, 0.0);
  const Vector2D e1(1.0, 0.0);
  const Vector2D e2(0.0, 1.0);
  const Vector2D vec(4.0, 3.0);
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm(e1), 1.0);
  EXPECT_DOUBLE_EQ(norm(e2), 1.0);
  EXPECT_DOUBLE_EQ(norm(vec), std::sqrt(4.0 * 4.0 + 3.0 * 3.0));
}

TEST(Vector2DTest, Norm2) {
  const Vector2D zero(0.0, 0.0);
  const Vector2D e1(1.0, 0.0);
  const Vector2D e2(0.0, 1.0);
  const Vector2D vec(4.0, 3.0);
  EXPECT_DOUBLE_EQ(norm2(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm2(e1), 1.0);
  EXPECT_DOUBLE_EQ(norm2(e2), 1.0);
  EXPECT_DOUBLE_EQ(norm2(vec), 4.0 * 4.0 + 3.0 * 3.0);
}

TEST(Vector2DTest, Abs) {
  const Vector2D vec1(4.0, 3.0);
  const Vector2D vec2(-4.0, 3.0);
  const Vector2D vec3(4.0, -3.0);
  const Vector2D vec4(-4.0, -3.0);
  const Vector2D abs1 = abs(vec1);
  EXPECT_DOUBLE_EQ(abs1.x(), 4.0);
  EXPECT_DOUBLE_EQ(abs1.y(), 3.0);
  const Vector2D abs2 = abs(vec2);
  EXPECT_DOUBLE_EQ(abs2.x(), 4.0);
  EXPECT_DOUBLE_EQ(abs2.y(), 3.0);
  const Vector2D abs3 = abs(vec3);
  EXPECT_DOUBLE_EQ(abs3.x(), 4.0);
  EXPECT_DOUBLE_EQ(abs3.y(), 3.0);
  const Vector2D abs4 = abs(vec4);
  EXPECT_DOUBLE_EQ(abs4.x(), 4.0);
  EXPECT_DOUBLE_EQ(abs4.y(), 3.0);
}

TEST(Vector2DTest, Exp) {
  const Vector2D vec(0.5, 1.0);
  const Vector2D result = exp(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::exp(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::exp(vec.y()));
}

TEST(Vector2DTest, Log) {
  const Vector2D vec(0.5, 1.0);
  const Vector2D result = log(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::log(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::log(vec.y()));
}

TEST(Vector2DTest, Sin) {
  const Vector2D vec(0, PI / 2);
  const Vector2D result = sin(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::sin(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::sin(vec.y()));
}

TEST(Vector2DTest, Cos) {
  const Vector2D vec(0, PI / 2);
  const Vector2D result = cos(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::cos(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::cos(vec.y()));
}

TEST(Vector2DTest, Tan) {
  const Vector2D vec(0, PI / 2);
  const Vector2D result = tan(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::tan(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::tan(vec.y()));
}

TEST(Vector2DTest, Asin) {
  const Vector2D vec(0, PI / 4);
  const Vector2D result = asin(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::asin(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::asin(vec.y()));
}

TEST(Vector2DTest, Acos) {
  const Vector2D vec(0, PI / 4);
  const Vector2D result = acos(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::acos(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::acos(vec.y()));
}

TEST(Vector2DTest, Atan) {
  const Vector2D vec(0, PI / 4);
  const Vector2D result = atan(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::atan(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::atan(vec.y()));
}

TEST(Vector2DTest, PowRealExponent) {
  const real_t exponent = 2.0;
  const Vector2D vec1(4.0, 3.0);
  const Vector2D vec2(-4.0, 3.0);
  const Vector2D result1 = pow(vec1, exponent);
  const Vector2D result2 = pow(vec1, exponent);
  EXPECT_DOUBLE_EQ(result1.x(), std::pow(vec1.x(), exponent));
  EXPECT_DOUBLE_EQ(result1.y(), std::pow(vec1.y(), exponent));
  EXPECT_DOUBLE_EQ(result2.x(), std::pow(vec2.x(), exponent));
  EXPECT_DOUBLE_EQ(result2.y(), std::pow(vec2.y(), exponent));
}

TEST(Vector2DTest, Sqrt) {
  const Vector2D vec1(1.0, 2.0);
  const Vector2D vec2(3.0, 4.0);
  const Vector2D vec3(3.0, -4.0);
  const Vector2D result1 = sqrt(vec1);
  const Vector2D result2 = sqrt(vec2);
  EXPECT_DOUBLE_EQ(result1.x(), std::sqrt(vec1.x()));
  EXPECT_DOUBLE_EQ(result1.y(), std::sqrt(vec1.y()));
  EXPECT_DOUBLE_EQ(result2.x(), std::sqrt(vec2.x()));
  EXPECT_DOUBLE_EQ(result2.y(), std::sqrt(vec2.y()));
  EXPECT_THROW(sqrt(vec3), mathlib::MathLibError);
}

TEST(Vector2DTest, Dot) {
  const Vector2D vec1(1, 2);
  const Vector2D vec2(3, 4);
  const Vector2D zero(0.0, 0.0);
  EXPECT_DOUBLE_EQ(dot(vec1, zero), 0.0);
  EXPECT_DOUBLE_EQ(dot(zero, vec1), 0.0);
  EXPECT_DOUBLE_EQ(dot(vec1, vec2), 1.0 * 3.0 + 2.0 * 4.0);
  EXPECT_DOUBLE_EQ(dot(vec2, vec1), 1.0 * 3.0 + 2.0 * 4.0);
}

TEST(Vector2DTest, Normalize) {
  const Vector2D zero(0.0, 0.0);
  EXPECT_THROW(normalize(zero), mathlib::DivByZeroError);

  const Vector2D e1(1.0, 0.0);
  const Vector2D n1 = normalize(e1);
  EXPECT_DOUBLE_EQ(n1.x(), 1.0);
  EXPECT_DOUBLE_EQ(n1.y(), 0.0);

  const Vector2D e2(0.0, 1.0);
  const Vector2D n2 = normalize(e2);
  EXPECT_DOUBLE_EQ(n2.x(), 0.0);
  EXPECT_DOUBLE_EQ(n2.y(), 1.0);

  const Vector2D vec(4.0, 2.0);
  const Vector2D n = normalize(vec);
  EXPECT_DOUBLE_EQ(n.x(), 4.0 / norm(vec));
  EXPECT_DOUBLE_EQ(n.y(), 2.0 / norm(vec));
}

TEST(Vector2DTest, Angle) {
  const Vector2D e1(1.0, 0.0);
  const Vector2D e2(0.0, 1.0);
  EXPECT_DOUBLE_EQ(angle(e1), 0.0);
  EXPECT_DOUBLE_EQ(angle(e2), PI/2);
  EXPECT_DOUBLE_EQ(angle(e1 + e2), PI/4);
}

TEST(Vector2DTest, AngleBetween2Vecs) {
  const Vector2D e1(1.0, 0.0);
  const Vector2D e2(0.0, 1.0);
  const Vector2D ones(1.0, 1.0);
  EXPECT_DOUBLE_EQ(angle(e1, e2), PI / 2);
  EXPECT_DOUBLE_EQ(angle(e1, ones), PI / 4);
  EXPECT_DOUBLE_EQ(angle(e1, -e1), -PI);
  EXPECT_DOUBLE_EQ(angle(-e2, e2), PI);
  EXPECT_DOUBLE_EQ(angle(e1, e1), 0);
}

TEST(Vector2DTest, AdditionVector2DVector2D) {
  const Vector2D vec1(1, 2);
  const Vector2D vec2(3, 4);
  const Vector2D result = vec1 + vec2;
  EXPECT_DOUBLE_EQ(result.x(), vec1.x() + vec2.x());
  EXPECT_DOUBLE_EQ(result.y(), vec1.y() + vec2.y());
}

TEST(Vector2DTest, AdditionVector2DReal) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = vec + value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() + value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() + value);
}

TEST(Vector2DTest, AdditionRealVector2D) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = value + vec;
  EXPECT_DOUBLE_EQ(result.x(), value + vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value + vec.y());
}

TEST(Vector2DTest, SubtractionVector2DVector2D) {
  const Vector2D vec1(1, 2);
  const Vector2D vec2(3, 4);
  const Vector2D result = vec1 - vec2;
  EXPECT_DOUBLE_EQ(result.x(), vec1.x() - vec2.x());
  EXPECT_DOUBLE_EQ(result.y(), vec1.y() - vec2.y());
}

TEST(Vector2DTest, SubtractionVector2DReal) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = vec - value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() - value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() - value);
}

TEST(Vector2DTest, SubtractionRealVector2D) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = value - vec;
  EXPECT_DOUBLE_EQ(result.x(), value - vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value - vec.y());
}

TEST(Vector2DTest, MultiplicationVector2DVector2D) {
  const Vector2D vec1(1, 2);
  const Vector2D vec2(3, 4);
  EXPECT_DOUBLE_EQ(vec1 * vec2, 1 * 3 + 2 * 4);
  EXPECT_DOUBLE_EQ(vec2 * vec1, 1 * 3 + 2 * 4);
}

TEST(Vector2DTest, MultiplicationVector2DReal) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = vec * value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() * value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() * value);
}

TEST(Vector2DTest, MultiplicationRealVector2D) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = value * vec;
  EXPECT_DOUBLE_EQ(result.x(), value * vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value * vec.y());
}

TEST(Vector2DTest, DivisionVector2DReal) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = vec / value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() / value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() / value);
  EXPECT_THROW(vec / 0.0, mathlib::DivByZeroError);
}

TEST(Vector2DTest, DivisionRealVector2D) {
  const real_t value = 2.0;
  const Vector2D vec(1, 2);
  const Vector2D result = value / vec;
  EXPECT_DOUBLE_EQ(result.x(), value / vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value / vec.y());
  EXPECT_THROW(value / Vector2D(0.0, 1.0), mathlib::DivByZeroError);
  EXPECT_THROW(value / Vector2D(1.0, 0.0), mathlib::DivByZeroError);
  EXPECT_THROW(value / Vector2D(0.0, 0.0), mathlib::DivByZeroError);
}