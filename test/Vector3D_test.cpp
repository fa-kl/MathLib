#include "DivByZeroError.hpp"
#include "MathLibError.hpp"
#include "Vector3D.hpp"

#include <complex>
#include <gtest/gtest.h>

using namespace mathlib;

TEST(Vector3DTest, DefaultConstructor) {
  const Vector3D vec;
  EXPECT_DOUBLE_EQ(vec.x(), 0.0);
  EXPECT_DOUBLE_EQ(vec.y(), 0.0);
  EXPECT_DOUBLE_EQ(vec.z(), 0.0);
  EXPECT_DOUBLE_EQ(x(vec), 0.0);
  EXPECT_DOUBLE_EQ(y(vec), 0.0);
  EXPECT_DOUBLE_EQ(z(vec), 0.0);
}

TEST(Vector3DTest, Constructor) {
  const Vector3D vec1(1.0, 2.0, 3.0);
  EXPECT_DOUBLE_EQ(vec1.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec1.y(), 2.0);
  EXPECT_DOUBLE_EQ(vec1.z(), 3.0);
  const Vector3D vec2(-2.0, -1.0, 5.0);
  EXPECT_DOUBLE_EQ(vec2.x(), -2.0);
  EXPECT_DOUBLE_EQ(vec2.y(), -1.0);
  EXPECT_DOUBLE_EQ(vec2.z(), 5.0);
}

TEST(Vector3DTest, CopyConstructor) {
  const Vector3D vec1(1.0, 2.0, 3.0);
  const Vector3D vec2(vec1);
  EXPECT_DOUBLE_EQ(vec2.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 2.0);
  EXPECT_DOUBLE_EQ(vec2.z(), 3.0);
}

TEST(Vector3DTest, Copy) {
  const Vector3D vec1(1.0, 2.0, 3.0);
  Vector3D vec2;
  vec2 = vec1;
  EXPECT_DOUBLE_EQ(vec2.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 2.0);
  EXPECT_DOUBLE_EQ(vec2.z(), 3.0);
}

TEST(Vector3DTest, Getter) {
  const Vector3D vec1(1.0, 2.0, 3.0);
  const Vector3D vec2(3.0, 4.0, 5.0);
  EXPECT_DOUBLE_EQ(vec1.x(), 1.0);
  EXPECT_DOUBLE_EQ(vec1.y(), 2.0);
  EXPECT_DOUBLE_EQ(vec1.z(), 3.0);
  EXPECT_DOUBLE_EQ(vec2.x(), 3.0);
  EXPECT_DOUBLE_EQ(vec2.y(), 4.0);
  EXPECT_DOUBLE_EQ(vec2.z(), 5.0);
}

TEST(Vector3DTest, Norm) {
  const Vector3D zero(0.0, 0.0, 0.0);
  const Vector3D e1(1.0, 0.0, 0.0);
  const Vector3D e2(0.0, 1.0, 0.0);
  const Vector3D e3(0.0, 0.0, 1.0);
  const Vector3D vec(4.0, 3.0, 12.0);
  EXPECT_DOUBLE_EQ(norm(zero), 0.0);
  EXPECT_DOUBLE_EQ(norm(e1), 1.0);
  EXPECT_DOUBLE_EQ(norm(e2), 1.0);
  EXPECT_DOUBLE_EQ(norm(e3), 1.0);
  EXPECT_DOUBLE_EQ(norm(vec), std::sqrt(4.0 * 4.0 + 3.0 * 3.0 + 12.0 * 12.0));
}

TEST(Vector3DTest, Norm2) {
  const Vector3D vec(4.0, 3.0, 12.0);
  EXPECT_DOUBLE_EQ(norm2(vec), 4.0 * 4.0 + 3.0 * 3.0 + 12.0 * 12.0);
}

TEST(Vector3DTest, Abs) {
  const Vector3D vec1(4.0, 3.0, -2.0);
  const Vector3D result = abs(vec1);
  EXPECT_DOUBLE_EQ(result.x(), 4.0);
  EXPECT_DOUBLE_EQ(result.y(), 3.0);
  EXPECT_DOUBLE_EQ(result.z(), 2.0);
}

TEST(Vector3DTest, Exp) {
  const Vector3D vec(0.5, 1.0, -1.0);
  const Vector3D result = exp(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::exp(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::exp(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::exp(vec.z()));
}

TEST(Vector3DTest, Log) {
  const Vector3D vec(0.5, 1.0, 2.0);
  const Vector3D result = log(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::log(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::log(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::log(vec.z()));
}

TEST(Vector3DTest, Sin) {
  const Vector3D vec(0, PI / 2, PI);
  const Vector3D result = sin(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::sin(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::sin(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::sin(vec.z()));
}

TEST(Vector3DTest, Cos) {
  const Vector3D vec(0, PI / 2, PI);
  const Vector3D result = cos(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::cos(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::cos(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::cos(vec.z()));
}

TEST(Vector3DTest, Tan) {
  const Vector3D vec(0, PI / 4, PI / 3);
  const Vector3D result = tan(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::tan(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::tan(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::tan(vec.z()));
}

TEST(Vector3DTest, Asin) {
  const Vector3D vec(0.0, 0.5, -0.5);
  const Vector3D result = asin(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::asin(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::asin(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::asin(vec.z()));
}

TEST(Vector3DTest, Acos) {
  const Vector3D vec(0.0, 0.5, -0.5);
  const Vector3D result = acos(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::acos(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::acos(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::acos(vec.z()));
}

TEST(Vector3DTest, Atan) {
  const Vector3D vec(0.0, 1.0, -1.0);
  const Vector3D result = atan(vec);
  EXPECT_DOUBLE_EQ(result.x(), std::atan(vec.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::atan(vec.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::atan(vec.z()));
}

TEST(Vector3DTest, PowRealExponent) {
  const real_t exponent = 2.0;
  const Vector3D vec(4.0, 3.0, 2.0);
  const Vector3D result = pow(vec, exponent);
  EXPECT_DOUBLE_EQ(result.x(), std::pow(vec.x(), exponent));
  EXPECT_DOUBLE_EQ(result.y(), std::pow(vec.y(), exponent));
  EXPECT_DOUBLE_EQ(result.z(), std::pow(vec.z(), exponent));
}

TEST(Vector3DTest, Sqrt) {
  const Vector3D vec1(1.0, 2.0, 3.0);
  const Vector3D result = sqrt(vec1);
  EXPECT_DOUBLE_EQ(result.x(), std::sqrt(vec1.x()));
  EXPECT_DOUBLE_EQ(result.y(), std::sqrt(vec1.y()));
  EXPECT_DOUBLE_EQ(result.z(), std::sqrt(vec1.z()));
  EXPECT_THROW(sqrt(Vector3D(-1.0, 2.0, 3.0)), mathlib::MathLibError);
}

TEST(Vector3DTest, Dot) {
  const Vector3D vec1(1, 2, 3);
  const Vector3D vec2(4, 5, 6);
  const Vector3D zero(0, 0, 0);
  EXPECT_DOUBLE_EQ(dot(vec1, zero), 0.0);
  EXPECT_DOUBLE_EQ(dot(zero, vec1), 0.0);
  EXPECT_DOUBLE_EQ(dot(vec1, vec2), 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0);
  EXPECT_DOUBLE_EQ(dot(vec2, vec1), 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0);
}

TEST(Vector3DTest, Normalize) {
  const Vector3D zero(0, 0, 0);
  EXPECT_THROW(normalize(zero), mathlib::DivByZeroError);

  const Vector3D e1(1.0, 0.0, 0.0);
  const Vector3D n1 = normalize(e1);
  EXPECT_DOUBLE_EQ(n1.x(), 1.0);
  EXPECT_DOUBLE_EQ(n1.y(), 0.0);
  EXPECT_DOUBLE_EQ(n1.z(), 0.0);

  const Vector3D e2(0.0, 1.0, 0.0);
  const Vector3D n2 = normalize(e2);
  EXPECT_DOUBLE_EQ(n2.x(), 0.0);
  EXPECT_DOUBLE_EQ(n2.y(), 1.0);
  EXPECT_DOUBLE_EQ(n2.z(), 0.0);

  const Vector3D vec(4.0, 2.0, -3.0);
  const Vector3D n = normalize(vec);
  EXPECT_DOUBLE_EQ(n.x(), 4.0 / norm(vec));
  EXPECT_DOUBLE_EQ(n.y(), 2.0 / norm(vec));
  EXPECT_DOUBLE_EQ(n.z(), -3.0 / norm(vec));
}

TEST(Vector3DTest, Cross) {
  const Vector3D e1(1.0, 0.0, 0.0);
  const Vector3D e2(0.0, 1.0, 0.0);
  const Vector3D e3(0.0, 0.0, 1.0);

  const Vector3D r1 = cross(e1, e2);
  EXPECT_DOUBLE_EQ(r1.x(), e3.x());
  EXPECT_DOUBLE_EQ(r1.y(), e3.y());
  EXPECT_DOUBLE_EQ(r1.z(), e3.z());

  const Vector3D r2 = cross(e1, e3);
  EXPECT_DOUBLE_EQ(r2.x(), -e2.x());
  EXPECT_DOUBLE_EQ(r2.y(), -e2.y());
  EXPECT_DOUBLE_EQ(r2.z(), -e2.z());

  const Vector3D r3 = cross(e2, e3);
  EXPECT_DOUBLE_EQ(r3.x(), e1.x());
  EXPECT_DOUBLE_EQ(r3.y(), e1.y());
  EXPECT_DOUBLE_EQ(r3.z(), e1.z());

  const Vector3D r4 = cross(e1, 2 * e1);
  EXPECT_DOUBLE_EQ(r4.x(), 0.0);
  EXPECT_DOUBLE_EQ(r4.y(), 0.0);
  EXPECT_DOUBLE_EQ(r4.z(), 0.0);
}

TEST(Vector3DTest, AdditionVector3DVector3D) {
  const Vector3D v1(1, 2, 3);
  const Vector3D v2(3, 4, 5);
  const Vector3D result = v1 + v2;
  EXPECT_DOUBLE_EQ(result.x(), v1.x() + v2.x());
  EXPECT_DOUBLE_EQ(result.y(), v1.y() + v2.y());
  EXPECT_DOUBLE_EQ(result.z(), v1.z() + v2.z());
}

TEST(Vector3DTest, AdditionVector3DReal) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = vec + value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() + value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() + value);
  EXPECT_DOUBLE_EQ(result.z(), vec.z() + value);
}

TEST(Vector3DTest, AdditionRealVector3D) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = value + vec;
  EXPECT_DOUBLE_EQ(result.x(), value + vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value + vec.y());
  EXPECT_DOUBLE_EQ(result.z(), value + vec.z());
}

TEST(Vector3DTest, SubtractionVector3DVector3D) {
  const Vector3D v1(1, 2, 3);
  const Vector3D v2(3, 4, 5);
  const Vector3D result = v1 - v2;
  EXPECT_DOUBLE_EQ(result.x(), v1.x() - v2.x());
  EXPECT_DOUBLE_EQ(result.y(), v1.y() - v2.y());
  EXPECT_DOUBLE_EQ(result.z(), v1.z() - v2.z());
}

TEST(Vector3DTest, SubtractionVector3DReal) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = vec - value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() - value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() - value);
  EXPECT_DOUBLE_EQ(result.z(), vec.z() - value);
}

TEST(Vector3DTest, SubtractionRealVector3D) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = value - vec;
  EXPECT_DOUBLE_EQ(result.x(), value - vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value - vec.y());
  EXPECT_DOUBLE_EQ(result.z(), value - vec.z());
}

TEST(Vector3DTest, MultiplicationVector3DVector3D) {
  const Vector3D v1(1, 2, 3);
  const Vector3D v2(3, 4, 5);
  EXPECT_DOUBLE_EQ(v1 * v2, 1 * 3 + 2 * 4 + 3 * 5);
  EXPECT_DOUBLE_EQ(v2 * v1, 1 * 3 + 2 * 4 + 3 * 5);
}

TEST(Vector3DTest, MultiplicationVector3DReal) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = vec * value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() * value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() * value);
  EXPECT_DOUBLE_EQ(result.z(), vec.z() * value);
}

TEST(Vector3DTest, MultiplicationRealVector3D) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = value * vec;
  EXPECT_DOUBLE_EQ(result.x(), value * vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value * vec.y());
  EXPECT_DOUBLE_EQ(result.z(), value * vec.z());
}

TEST(Vector3DTest, DivisionVector3DReal) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = vec / value;
  EXPECT_DOUBLE_EQ(result.x(), vec.x() / value);
  EXPECT_DOUBLE_EQ(result.y(), vec.y() / value);
  EXPECT_DOUBLE_EQ(result.z(), vec.z() / value);
  EXPECT_THROW(vec / 0.0, mathlib::DivByZeroError);
}

TEST(Vector3DTest, DivisionRealVector3D) {
  const real_t value = 2.0;
  const Vector3D vec(1, 2, 3);
  const Vector3D result = value / vec;
  EXPECT_DOUBLE_EQ(result.x(), value / vec.x());
  EXPECT_DOUBLE_EQ(result.y(), value / vec.y());
  EXPECT_DOUBLE_EQ(result.z(), value / vec.z());
  EXPECT_THROW(value / Vector3D(0.0, 1.0, 1.0), mathlib::DivByZeroError);
  EXPECT_THROW(value / Vector3D(1.0, 0.0, 1.0), mathlib::DivByZeroError);
  EXPECT_THROW(value / Vector3D(1.0, 1.0, 0.0), mathlib::DivByZeroError);
}
