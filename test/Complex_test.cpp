#include "Complex.hpp"

#include <complex>
#include <gtest/gtest.h>

using namespace mathlib;

constexpr real_t TEST_EPSILON = 1e-9;

TEST(ComplexTest, DefaultConstructor) {
  const Complex z;
  const std::complex<real_t> std_z;
  EXPECT_DOUBLE_EQ(z.real(), std_z.real());
  EXPECT_DOUBLE_EQ(z.imag(), std_z.imag());
}

TEST(ComplexTest, Constructor) {
  const Complex z1(1, 2);
  const std::complex<real_t> std_z1(1, 2);
  EXPECT_DOUBLE_EQ(z1.real(), std_z1.real());
  EXPECT_DOUBLE_EQ(z1.imag(), std_z1.imag());

  const Complex z2(-1, -2);
  const std::complex<real_t> std_z2(-1, -2);
  EXPECT_DOUBLE_EQ(z2.real(), std_z2.real());
  EXPECT_DOUBLE_EQ(z2.imag(), std_z2.imag());
}

TEST(ComplexTest, CopyConstructor) {
  const Complex z1(1, 2);
  const Complex z2(z1);
  const std::complex<real_t> std_z1(1, 2);
  const std::complex<real_t> std_z2(std_z1);
  EXPECT_DOUBLE_EQ(z2.real(), std_z2.real());
  EXPECT_DOUBLE_EQ(z2.imag(), std_z2.imag());
}

TEST(ComplexTest, Copy) {
  const Complex z1(2, 1);
  Complex z2;
  z2 = z1;
  const std::complex<real_t> std_z1(2, 1);
  const std::complex<real_t> std_z2 = std_z1;
  EXPECT_DOUBLE_EQ(z2.real(), std_z2.real());
  EXPECT_DOUBLE_EQ(z2.imag(), std_z2.imag());
}

TEST(ComplexTest, Real) {
  const Complex z(3.5, -2.1);
  const std::complex<real_t> std_z(3.5, -2.1);
  EXPECT_DOUBLE_EQ(z.real(), std_z.real());
  EXPECT_DOUBLE_EQ(real(z), std_z.real());
}

TEST(ComplexTest, Imag) {
  const Complex z(3.5, -2.1);
  const std::complex<real_t> std_z(3.5, -2.1);
  EXPECT_DOUBLE_EQ(z.imag(), std_z.imag());
  EXPECT_DOUBLE_EQ(imag(z), std_z.imag());
}

TEST(ComplexTest, Abs) {
  const Complex z1(3, 4);
  const std::complex<real_t> std_z1(3, 4);
  EXPECT_DOUBLE_EQ(abs(z1), std::abs(std_z1));

  const Complex z2(-1, -1);
  const std::complex<real_t> std_z2(-1, -1);
  EXPECT_DOUBLE_EQ(abs(z2), std::abs(std_z2));
}

TEST(ComplexTest, Abs2) {
  const Complex z1(3, 4);
  const std::complex<real_t> std_z1(3, 4);
  EXPECT_DOUBLE_EQ(abs2(z1), std::norm(std_z1));

  const Complex z2(-1, -1);
  const std::complex<real_t> std_z2(-1, -1);
  EXPECT_DOUBLE_EQ(abs2(z2), std::norm(std_z2));
}

TEST(ComplexTest, Arg) {
  const Complex z1(1, 1);
  const std::complex<real_t> std_z1(1, 1);
  EXPECT_DOUBLE_EQ(arg(z1), std::arg(std_z1));

  const Complex z2(-1, 0);
  const std::complex<real_t> std_z2(-1, 0);
  EXPECT_DOUBLE_EQ(arg(z2), std::arg(std_z2));
}

TEST(ComplexTest, Conj) {
  const Complex z(2, -3);
  const std::complex<real_t> std_z(2, -3);
  const Complex result = conj(z);
  const std::complex<real_t> std_result = std::conj(std_z);
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, Exp) {
  const Complex z(1, 0.5);
  const std::complex<real_t> std_z(1, 0.5);
  const Complex result = exp(z);
  const std::complex<real_t> std_result = std::exp(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Log) {
  const Complex z(2, 1);
  const std::complex<real_t> std_z(2, 1);
  const Complex result = log(z);
  const std::complex<real_t> std_result = std::log(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Sin) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = sin(z);
  const std::complex<real_t> std_result = std::sin(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Cos) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = cos(z);
  const std::complex<real_t> std_result = std::cos(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Tan) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = tan(z);
  const std::complex<real_t> std_result = std::tan(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Asin) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = asin(z);
  const std::complex<real_t> std_result = std::asin(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Acos) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = acos(z);
  const std::complex<real_t> std_result = std::acos(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Atan) {
  const Complex z(0.5, 0.3);
  const std::complex<real_t> std_z(0.5, 0.3);
  const Complex result = atan(z);
  const std::complex<real_t> std_result = std::atan(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, PowRealExponent) {
  const Complex z(4, 3);
  const real_t w = 2.5;
  const std::complex<real_t> std_z(4, 3);
  const std::complex<real_t> std_result = std::pow(std_z, w);
  const Complex result = pow(z, w);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, PowComplexExponent) {
  const Complex z(4, 3);
  const Complex w(1.2, -0.7);
  const std::complex<real_t> std_z(4, 3);
  const std::complex<real_t> std_w(1.2, -0.7);
  const std::complex<real_t> std_result = std::pow(std_z, std_w);
  const Complex result = pow(z, w);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, Sqrt) {
  const Complex z(4, 3);
  const std::complex<real_t> std_z(4, 3);
  const Complex result = sqrt(z);
  const std::complex<real_t> std_result = std::sqrt(std_z);
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, AdditionComplexComplex) {
  const Complex z1(1, 2);
  const Complex z2(3, 4);
  const std::complex<real_t> std_z1(1, 2);
  const std::complex<real_t> std_z2(3, 4);
  const Complex result = z1 + z2;
  const std::complex<real_t> std_result = std_z1 + std_z2;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, AdditionComplexReal) {
  const Complex z(1, 2);
  const real_t r = 3;
  const std::complex<real_t> std_z(1, 2);
  const Complex result = z + r;
  const std::complex<real_t> std_result = std_z + r;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, AdditionRealComplex) {
  const real_t r = 3;
  const Complex z(1, 2);
  const std::complex<real_t> std_z(1, 2);
  const Complex result = r + z;
  const std::complex<real_t> std_result = r + std_z;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, SubtractionComplexComplex) {
  const Complex z1(5, 6);
  const Complex z2(1, 2);
  const std::complex<real_t> std_z1(5, 6);
  const std::complex<real_t> std_z2(1, 2);
  const Complex result = z1 - z2;
  const std::complex<real_t> std_result = std_z1 - std_z2;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, SubtractionComplexReal) {
  const Complex z(5, 2);
  const real_t r = 3;
  const std::complex<real_t> std_z(5, 2);
  const Complex result = z - r;
  const std::complex<real_t> std_result = std_z - r;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, SubtractionRealComplex) {
  const real_t r = 5;
  const Complex z(1, 2);
  const std::complex<real_t> std_z(1, 2);
  const Complex result = r - z;
  const std::complex<real_t> std_result = r - std_z;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, MultiplicationComplexComplex) {
  const Complex z1(1, 2);
  const Complex z2(3, 4);
  const std::complex<real_t> std_z1(1, 2);
  const std::complex<real_t> std_z2(3, 4);
  const Complex result = z1 * z2;
  const std::complex<real_t> std_result = std_z1 * std_z2;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, MultiplicationComplexReal) {
  const Complex z(1, 2);
  const real_t r = 3;
  const std::complex<real_t> std_z(1, 2);
  const Complex result = z * r;
  const std::complex<real_t> std_result = std_z * r;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, MultiplicationRealComplex) {
  const real_t r = 3;
  const Complex z(1, 2);
  const std::complex<real_t> std_z(1, 2);
  const Complex result = r * z;
  const std::complex<real_t> std_result = r * std_z;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, DivisionComplexComplex) {
  const Complex z1(10, 5);
  const Complex z2(2, 1);
  const std::complex<real_t> std_z1(10, 5);
  const std::complex<real_t> std_z2(2, 1);
  const Complex result = z1 / z2;
  const std::complex<real_t> std_result = std_z1 / std_z2;
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}

TEST(ComplexTest, DivisionComplexReal) {
  const Complex z(6, 4);
  const real_t r = 2;
  const std::complex<real_t> std_z(6, 4);
  const Complex result = z / r;
  const std::complex<real_t> std_result = std_z / r;
  EXPECT_DOUBLE_EQ(result.real(), std_result.real());
  EXPECT_DOUBLE_EQ(result.imag(), std_result.imag());
}

TEST(ComplexTest, DivisionRealComplex) {
  const real_t r = 6;
  const Complex z(2, 1);
  const std::complex<real_t> std_z(2, 1);
  const Complex result = r / z;
  const std::complex<real_t> std_result = r / std_z;
  EXPECT_NEAR(result.real(), std_result.real(), TEST_EPSILON);
  EXPECT_NEAR(result.imag(), std_result.imag(), TEST_EPSILON);
}