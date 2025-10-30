#include "rounding.hpp"
#include <gtest/gtest.h>
#include <iostream>

TEST(RoundingTest, IntegerRounding) {
  EXPECT_TRUE(mathlib::round(3.4f) == 3);
  EXPECT_TRUE(mathlib::round(3.5f) == 4);
  EXPECT_TRUE(mathlib::round(-3.5f) == -4);
  EXPECT_TRUE(mathlib::round(2.0f, mathlib::CEIL) == 2);
  EXPECT_TRUE(mathlib::round(2.1f, mathlib::CEIL) == 3);
  EXPECT_TRUE(mathlib::round(2.9f, mathlib::FLOOR) == 2);
  EXPECT_TRUE(mathlib::round(-2.1f, mathlib::CEIL) == -2);
  EXPECT_TRUE(mathlib::round(-2.1f, mathlib::FLOOR) == -3);
}

TEST(RoundingTest, NearestRounding) {
  EXPECT_TRUE(abs(mathlib::nearest(3.14159f, 2) - 3.14f) < 1e-6);
  EXPECT_TRUE(abs(mathlib::nearest(3.145f, 2) - 3.15f) < 1e-6);
  EXPECT_TRUE(abs(mathlib::nearest(-3.145f, 2) + 3.15f) < 1e-6);
}

TEST(RoundingTest, CeilingRounding) {
  EXPECT_TRUE(abs(mathlib::ceil(3.141f, 2) - 3.15f) < 1e-6);
  EXPECT_TRUE(abs(mathlib::ceil(-3.141f, 2) + 3.14f) < 1e-6);
}

TEST(RoundingTest, FloorRounding) {
  EXPECT_TRUE(abs(mathlib::floor(3.149f, 2) - 3.14f) < 1e-6);
  EXPECT_TRUE(abs(mathlib::floor(-3.149f, 2) + 3.15f) < 1e-6);
}

TEST(RoundingTest, DoubleRounding) {
  EXPECT_TRUE(abs(mathlib::nearest(3.14159265358979, 3) - 3.142) < 1e-12);
  EXPECT_TRUE(abs(mathlib::ceil(3.14159265358979, 3) - 3.142) < 1e-12);
  EXPECT_TRUE(abs(mathlib::floor(3.14159265358979, 3) - 3.141) < 1e-12);
}

TEST(RoundingTest, ExtremeDecimals) {
  EXPECT_TRUE(abs(mathlib::nearest(1.23456789, 6) - 1.234568) < 1e-9);
  EXPECT_TRUE(abs(mathlib::floor(1.23456789, 6) - 1.234567) < 1e-9);
}

TEST(RoundingTest, LargeScaleRounding) {
  EXPECT_TRUE(abs(mathlib::nearest(123456.789, 0) - 123457) < 1e-6);
  EXPECT_TRUE(abs(mathlib::floor(123456.789, 0) - 123456) < 1e-6);
  EXPECT_TRUE(abs(mathlib::ceil(123456.001, 0) - 123457) < 1e-6);
}

TEST(RoundingTest, NegativeRoundingChecks) {
  EXPECT_TRUE(abs(mathlib::nearest(-1.2345, 2) + 1.23) < 1e-6);
  EXPECT_TRUE(abs(mathlib::ceil(-1.2345, 2) + 1.23) < 1e-6);
  EXPECT_TRUE(abs(mathlib::floor(-1.2345, 2) + 1.24) < 1e-6);
}

TEST(RoundingTest, ZeroHandling) {
  EXPECT_TRUE(mathlib::nearest(0.0f, 3) == 0.0f);
  EXPECT_TRUE(mathlib::ceil(0.0f, 3) == 0.0f);
  EXPECT_TRUE(mathlib::floor(0.0f, 3) == 0.0f);
}
