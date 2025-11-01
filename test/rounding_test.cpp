/*****************************************************************************************
 * @file: rounding_test.cpp
 *
 * @brief: Unit tests for rounding functionality
 *
 * @author: fakl
 * @date: November 2025
 *
 ****************************************************************************************/

#include "rounding.hpp"

#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "config.hpp"

using namespace mathlib;

// ======================================================================================
// Integer Rounding Tests (no decimal precision)
// ======================================================================================

TEST(RoundingTest, IntegerRounding_Nearest)
{
  // Positive values
  EXPECT_EQ(round(1.4f, NEAREST), 1);
  EXPECT_EQ(round(1.5f, NEAREST), 2);
  EXPECT_EQ(round(1.6f, NEAREST), 2);
  EXPECT_EQ(round(2.5f, NEAREST), 3);  // round half away from zero

  // Negative values
  EXPECT_EQ(round(-1.4f, NEAREST), -1);
  EXPECT_EQ(round(-1.5f, NEAREST), -2);
  EXPECT_EQ(round(-1.6f, NEAREST), -2);

  // Edge cases
  EXPECT_EQ(round(0.0f, NEAREST), 0);
  EXPECT_EQ(round(0.5f, NEAREST), 1);    // round half away from zero
  EXPECT_EQ(round(-0.5f, NEAREST), -1);  // round half away from zero
}

TEST(RoundingTest, IntegerRounding_Ceil)
{
  // Positive values
  EXPECT_EQ(round(1.1f, CEIL), 2);
  EXPECT_EQ(round(1.5f, CEIL), 2);
  EXPECT_EQ(round(1.9f, CEIL), 2);
  EXPECT_EQ(round(2.0f, CEIL), 2);

  // Negative values
  EXPECT_EQ(round(-1.1f, CEIL), -1);
  EXPECT_EQ(round(-1.5f, CEIL), -1);
  EXPECT_EQ(round(-1.9f, CEIL), -1);
  EXPECT_EQ(round(-2.0f, CEIL), -2);

  // Edge cases
  EXPECT_EQ(round(0.0f, CEIL), 0);
  EXPECT_EQ(round(0.1f, CEIL), 1);
  EXPECT_EQ(round(-0.1f, CEIL), 0);
}

TEST(RoundingTest, IntegerRounding_Floor)
{
  // Positive values
  EXPECT_EQ(round(1.1f, FLOOR), 1);
  EXPECT_EQ(round(1.5f, FLOOR), 1);
  EXPECT_EQ(round(1.9f, FLOOR), 1);
  EXPECT_EQ(round(2.0f, FLOOR), 2);

  // Negative values
  EXPECT_EQ(round(-1.1f, FLOOR), -2);
  EXPECT_EQ(round(-1.5f, FLOOR), -2);
  EXPECT_EQ(round(-1.9f, FLOOR), -2);
  EXPECT_EQ(round(-2.0f, FLOOR), -2);

  // Edge cases
  EXPECT_EQ(round(0.0f, FLOOR), 0);
  EXPECT_EQ(round(0.1f, FLOOR), 0);
  EXPECT_EQ(round(-0.1f, FLOOR), -1);
}

// ======================================================================================
// Float (real32_t) Rounding with Decimal Precision
// ======================================================================================

TEST(RoundingTest, Float_Nearest_Precision)
{
  real32_t value = 1.23456f;

  // Various precisions
  EXPECT_FLOAT_EQ(round(value, 0, NEAREST), 1.0f);
  EXPECT_FLOAT_EQ(round(value, 1, NEAREST), 1.2f);
  EXPECT_FLOAT_EQ(round(value, 2, NEAREST), 1.23f);
  EXPECT_FLOAT_EQ(round(value, 3, NEAREST), 1.235f);
  EXPECT_FLOAT_EQ(round(value, 4, NEAREST), 1.2346f);
  EXPECT_FLOAT_EQ(round(value, 5, NEAREST), 1.23456f);

  // Negative values
  value = -3.87654f;
  EXPECT_FLOAT_EQ(round(value, 0, NEAREST), -4.0f);
  EXPECT_FLOAT_EQ(round(value, 1, NEAREST), -3.9f);
  EXPECT_FLOAT_EQ(round(value, 2, NEAREST), -3.88f);
  EXPECT_FLOAT_EQ(round(value, 3, NEAREST), -3.877f);
}

TEST(RoundingTest, Float_Ceil_Precision)
{
  real32_t value = 1.23456f;

  EXPECT_FLOAT_EQ(round(value, 0, CEIL), 2.0f);
  EXPECT_FLOAT_EQ(round(value, 1, CEIL), 1.3f);
  EXPECT_FLOAT_EQ(round(value, 2, CEIL), 1.24f);
  EXPECT_FLOAT_EQ(round(value, 3, CEIL), 1.235f);

  // Negative values
  value = -1.23456f;
  EXPECT_FLOAT_EQ(round(value, 0, CEIL), -1.0f);
  EXPECT_FLOAT_EQ(round(value, 1, CEIL), -1.2f);
  EXPECT_FLOAT_EQ(round(value, 2, CEIL), -1.23f);
}

TEST(RoundingTest, Float_Floor_Precision)
{
  real32_t value = 1.23456f;

  EXPECT_FLOAT_EQ(round(value, 0, FLOOR), 1.0f);
  EXPECT_FLOAT_EQ(round(value, 1, FLOOR), 1.2f);
  EXPECT_FLOAT_EQ(round(value, 2, FLOOR), 1.23f);
  EXPECT_FLOAT_EQ(round(value, 3, FLOOR), 1.234f);

  // Negative values
  value = -1.23456f;
  EXPECT_FLOAT_EQ(round(value, 0, FLOOR), -2.0f);
  EXPECT_FLOAT_EQ(round(value, 1, FLOOR), -1.3f);
  EXPECT_FLOAT_EQ(round(value, 2, FLOOR), -1.24f);
}

TEST(RoundingTest, Float_EdgeCases)
{
  // Zero
  EXPECT_FLOAT_EQ(round(0.0f, 2, NEAREST), 0.0f);
  EXPECT_FLOAT_EQ(round(0.0f, 2, CEIL), 0.0f);
  EXPECT_FLOAT_EQ(round(0.0f, 2, FLOOR), 0.0f);

  // Very small values
  EXPECT_FLOAT_EQ(round(0.001f, 2, NEAREST), 0.0f);
  EXPECT_FLOAT_EQ(round(0.001f, 2, CEIL), 0.01f);
  EXPECT_FLOAT_EQ(round(0.001f, 2, FLOOR), 0.0f);

  // Exact values (no rounding needed)
  EXPECT_FLOAT_EQ(round(1.5f, 1, NEAREST), 1.5f);
  EXPECT_FLOAT_EQ(round(2.25f, 2, NEAREST), 2.25f);
}

// ======================================================================================
// Double (real64_t) Rounding with Decimal Precision
// ======================================================================================

TEST(RoundingTest, Double_Nearest_Precision)
{
  real64_t value = 1.23456789;

  // Various precisions
  EXPECT_DOUBLE_EQ(round(value, 0, NEAREST), 1.0);
  EXPECT_DOUBLE_EQ(round(value, 1, NEAREST), 1.2);
  EXPECT_DOUBLE_EQ(round(value, 2, NEAREST), 1.23);
  EXPECT_DOUBLE_EQ(round(value, 3, NEAREST), 1.235);
  EXPECT_DOUBLE_EQ(round(value, 4, NEAREST), 1.2346);
  EXPECT_DOUBLE_EQ(round(value, 5, NEAREST), 1.23457);
  EXPECT_DOUBLE_EQ(round(value, 6, NEAREST), 1.234568);

  // Negative values
  value = -9.87654321;
  EXPECT_DOUBLE_EQ(round(value, 0, NEAREST), -10.0);
  EXPECT_DOUBLE_EQ(round(value, 1, NEAREST), -9.9);
  EXPECT_DOUBLE_EQ(round(value, 2, NEAREST), -9.88);
  EXPECT_DOUBLE_EQ(round(value, 3, NEAREST), -9.877);
}

TEST(RoundingTest, Double_Ceil_Precision)
{
  real64_t value = 1.23456789;

  EXPECT_DOUBLE_EQ(round(value, 0, CEIL), 2.0);
  EXPECT_DOUBLE_EQ(round(value, 1, CEIL), 1.3);
  EXPECT_DOUBLE_EQ(round(value, 2, CEIL), 1.24);
  EXPECT_DOUBLE_EQ(round(value, 3, CEIL), 1.235);

  // Negative values
  value = -1.23456789;
  EXPECT_DOUBLE_EQ(round(value, 0, CEIL), -1.0);
  EXPECT_DOUBLE_EQ(round(value, 1, CEIL), -1.2);
  EXPECT_DOUBLE_EQ(round(value, 2, CEIL), -1.23);
}

TEST(RoundingTest, Double_Floor_Precision)
{
  real64_t value = 1.23456789;

  EXPECT_DOUBLE_EQ(round(value, 0, FLOOR), 1.0);
  EXPECT_DOUBLE_EQ(round(value, 1, FLOOR), 1.2);
  EXPECT_DOUBLE_EQ(round(value, 2, FLOOR), 1.23);
  EXPECT_DOUBLE_EQ(round(value, 3, FLOOR), 1.234);

  // Negative values
  value = -1.23456789;
  EXPECT_DOUBLE_EQ(round(value, 0, FLOOR), -2.0);
  EXPECT_DOUBLE_EQ(round(value, 1, FLOOR), -1.3);
  EXPECT_DOUBLE_EQ(round(value, 2, FLOOR), -1.24);
}

TEST(RoundingTest, Double_EdgeCases)
{
  // Zero
  EXPECT_DOUBLE_EQ(round(0.0, 2, NEAREST), 0.0);
  EXPECT_DOUBLE_EQ(round(0.0, 2, CEIL), 0.0);
  EXPECT_DOUBLE_EQ(round(0.0, 2, FLOOR), 0.0);

  // Very small values
  EXPECT_DOUBLE_EQ(round(0.001, 2, NEAREST), 0.0);
  EXPECT_DOUBLE_EQ(round(0.001, 2, CEIL), 0.01);
  EXPECT_DOUBLE_EQ(round(0.001, 2, FLOOR), 0.0);

  // Exact values (no rounding needed)
  EXPECT_DOUBLE_EQ(round(1.5, 1, NEAREST), 1.5);
  EXPECT_DOUBLE_EQ(round(2.25, 2, NEAREST), 2.25);

  // Large values
  EXPECT_DOUBLE_EQ(round(123456.789, 2, NEAREST), 123456.79);
  EXPECT_DOUBLE_EQ(round(123456.789, 1, CEIL), 123456.8);
  EXPECT_DOUBLE_EQ(round(123456.789, 1, FLOOR), 123456.7);
}

TEST(RoundingTest, Double_HighPrecision)
{
  real64_t value = 1.23456789012345;

  EXPECT_DOUBLE_EQ(round(value, 7, NEAREST), 1.2345679);
  EXPECT_DOUBLE_EQ(round(value, 8, NEAREST), 1.23456789);
  EXPECT_DOUBLE_EQ(round(value, 9, NEAREST), 1.234567890);
  EXPECT_DOUBLE_EQ(round(value, 10, NEAREST), 1.2345678901);
}

// ======================================================================================
// Convenience Wrapper Tests (Float)
// ======================================================================================

TEST(RoundingTest, Float_NearestWrapper)
{
  real32_t value = 1.23456f;

  EXPECT_FLOAT_EQ(nearest(value, 2), round(value, 2, NEAREST));
  EXPECT_FLOAT_EQ(nearest(value, 3), round(value, 3, NEAREST));
  EXPECT_FLOAT_EQ(nearest(-5.6789f, 2), round(-5.6789f, 2, NEAREST));
}

TEST(RoundingTest, Float_CeilWrapper)
{
  real32_t value = 1.23456f;

  EXPECT_FLOAT_EQ(mathlib::ceil(value, 2), round(value, 2, CEIL));
  EXPECT_FLOAT_EQ(mathlib::ceil(value, 3), round(value, 3, CEIL));
  EXPECT_FLOAT_EQ(mathlib::ceil(-5.6789f, 2), round(-5.6789f, 2, CEIL));
}

TEST(RoundingTest, Float_FloorWrapper)
{
  real32_t value = 1.23456f;

  EXPECT_FLOAT_EQ(mathlib::floor(value, 2), round(value, 2, FLOOR));
  EXPECT_FLOAT_EQ(mathlib::floor(value, 3), round(value, 3, FLOOR));
  EXPECT_FLOAT_EQ(mathlib::floor(-5.6789f, 2), round(-5.6789f, 2, FLOOR));
}

// ======================================================================================
// Convenience Wrapper Tests (Double)
// ======================================================================================

TEST(RoundingTest, Double_NearestWrapper)
{
  real64_t value = 1.23456789;

  EXPECT_DOUBLE_EQ(nearest(value, 2), round(value, 2, NEAREST));
  EXPECT_DOUBLE_EQ(nearest(value, 4), round(value, 4, NEAREST));
  EXPECT_DOUBLE_EQ(nearest(-5.6789012, 3), round(-5.6789012, 3, NEAREST));
}

TEST(RoundingTest, Double_CeilWrapper)
{
  real64_t value = 1.23456789;

  EXPECT_DOUBLE_EQ(mathlib::ceil(value, 2), round(value, 2, CEIL));
  EXPECT_DOUBLE_EQ(mathlib::ceil(value, 4), round(value, 4, CEIL));
  EXPECT_DOUBLE_EQ(mathlib::ceil(-5.6789012, 3), round(-5.6789012, 3, CEIL));
}

TEST(RoundingTest, Double_FloorWrapper)
{
  real64_t value = 1.23456789;

  EXPECT_DOUBLE_EQ(mathlib::floor(value, 2), round(value, 2, FLOOR));
  EXPECT_DOUBLE_EQ(mathlib::floor(value, 4), round(value, 4, FLOOR));
  EXPECT_DOUBLE_EQ(mathlib::floor(-5.6789012, 3), round(-5.6789012, 3, FLOOR));
}

// ======================================================================================
// Special Cases and Boundary Tests
// ======================================================================================

TEST(RoundingTest, ZeroPrecision)
{
  // Zero precision should round to nearest integer
  EXPECT_FLOAT_EQ(round(1.23f, 0, NEAREST), 1.0f);
  EXPECT_FLOAT_EQ(round(1.67f, 0, NEAREST), 2.0f);

  EXPECT_DOUBLE_EQ(round(9.87, 0, NEAREST), 10.0);
  EXPECT_DOUBLE_EQ(round(5.23, 0, NEAREST), 5.0);
}

TEST(RoundingTest, NegativeZeroHandling)
{
  // Ensure proper handling around zero
  EXPECT_FLOAT_EQ(round(-0.001f, 2, NEAREST), 0.0f);
  EXPECT_FLOAT_EQ(round(0.001f, 2, NEAREST), 0.0f);

  EXPECT_DOUBLE_EQ(round(-0.001, 2, NEAREST), 0.0);
  EXPECT_DOUBLE_EQ(round(0.001, 2, NEAREST), 0.0);
}

TEST(RoundingTest, RoundingHalfValues)
{
  // Test .5 values (round half away from zero)
  EXPECT_FLOAT_EQ(round(0.5f, 0, NEAREST), 1.0f);
  EXPECT_FLOAT_EQ(round(1.5f, 0, NEAREST), 2.0f);
  EXPECT_FLOAT_EQ(round(2.5f, 0, NEAREST), 3.0f);
  EXPECT_FLOAT_EQ(round(3.5f, 0, NEAREST), 4.0f);

  // With precision
  EXPECT_DOUBLE_EQ(round(1.25, 1, NEAREST), 1.3);
  EXPECT_DOUBLE_EQ(round(1.35, 1, NEAREST), 1.4);
}
