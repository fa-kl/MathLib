#include "compare.hpp"

#include <iostream>

#include <gtest/gtest.h>

using namespace mathlib;

TEST(CompareTest, isFuzzyEqualIdentical)
{
  EXPECT_TRUE(isFuzzyEqual(5.0, 5.0, 0.001));
  EXPECT_TRUE(isFuzzyEqual(0.0, 0.0, 0.001));
  EXPECT_TRUE(isFuzzyEqual(-3.0, -3.0, 0.001));
}

TEST(CompareTest, isFuzzyEqualWithinEpsilon)
{
  real_t epsilon = 0.001;
  EXPECT_TRUE(isFuzzyEqual(5.0 + 0.9 * epsilon, 5.0, epsilon));
  EXPECT_TRUE(isFuzzyEqual(5.0, 5.0 - 0.9 * epsilon, epsilon));
  EXPECT_TRUE(isFuzzyEqual(5.0 + 0.5 * epsilon, 5.0, epsilon));
}

TEST(CompareTest, isFuzzyEqualBeyondEpsilon)
{
  EXPECT_FALSE(isFuzzyEqual(5.0, 5.00101, 0.001));
  EXPECT_FALSE(isFuzzyEqual(5.0, 4.99899, 0.001));
  EXPECT_FALSE(isFuzzyEqual(5.0, 5.1, 0.001));
  EXPECT_FALSE(isFuzzyEqual(5.0, 4.9, 0.001));
}

TEST(CompareTest, isFuzzyGreaterIdentical)
{
  EXPECT_FALSE(isFuzzyGreater(5.0, 5.0, 0.001));
  EXPECT_FALSE(isFuzzyGreater(0.0, 0.0, 0.001));
}

TEST(CompareTest, isFuzzyGreaterWithinEpsilon)
{
  EXPECT_TRUE(isFuzzyGreater(5.001, 5.0, 0.001));
  EXPECT_TRUE(isFuzzyGreater(5.0, 4.999, 0.001));
  EXPECT_FALSE(isFuzzyGreater(5.0005, 5.0, 0.001));
}

TEST(CompareTest, isFuzzyGreaterBeyondEpsilon)
{
  EXPECT_TRUE(isFuzzyGreater(5.00101, 5.0, 0.001));
  EXPECT_TRUE(isFuzzyGreater(5.0, 4.99899, 0.001));
  EXPECT_TRUE(isFuzzyGreater(6.0, 5.0, 0.001));
}

TEST(CompareTest, isFuzzySmallerIdentical)
{
  EXPECT_FALSE(isFuzzySmaller(5.0, 5.0, 0.001));
  EXPECT_FALSE(isFuzzySmaller(0.0, 0.0, 0.001));
}

TEST(CompareTest, isFuzzySmallerWithinEpsilon)
{
  EXPECT_TRUE(isFuzzySmaller(4.999, 5.0, 0.001));
  EXPECT_TRUE(isFuzzySmaller(5.0, 5.001, 0.001));
  EXPECT_FALSE(isFuzzySmaller(4.9995, 5.0, 0.001));
}

TEST(CompareTest, isFuzzySmallerBeyondEpsilon)
{
  EXPECT_TRUE(isFuzzySmaller(4.99899, 5.0, 0.001));
  EXPECT_TRUE(isFuzzySmaller(5.0, 5.00101, 0.001));
  EXPECT_TRUE(isFuzzySmaller(4.0, 5.0, 0.001));
}

TEST(CompareTest, isStrictFuzzyGreaterIdentical)
{
  EXPECT_FALSE(isStrictFuzzyGreater(5.0, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzyGreater(0.0, 0.0, 0.001));
}

TEST(CompareTest, isStrictFuzzyGreaterWithinEpsilon)
{
  EXPECT_FALSE(isStrictFuzzyGreater(5.001, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzyGreater(5.0, 4.999, 0.001));
  EXPECT_FALSE(isStrictFuzzyGreater(5.0005, 5.0, 0.001));
}

TEST(CompareTest, isStrictFuzzyGreaterBeyondEpsilon)
{
  EXPECT_TRUE(isStrictFuzzyGreater(5.0 + 0.00101, 5.0, 0.001));
  EXPECT_TRUE(isStrictFuzzyGreater(6.0, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzyGreater(5.0, 6.0, 0.001));
  EXPECT_FALSE(isStrictFuzzyGreater(5.0 - 0.00101, 5.0, 0.001));
}

TEST(CompareTest, isStrictFuzzySmallerIdentical)
{
  EXPECT_FALSE(isStrictFuzzySmaller(5.0, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzySmaller(0.0, 0.0, 0.001));
}

TEST(CompareTest, isStrictFuzzySmallerWithinEpsilon)
{
  EXPECT_FALSE(isStrictFuzzySmaller(5.0 - 0.001, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzySmaller(5.0, 5.0 + 0.001, 0.001));
  EXPECT_FALSE(isStrictFuzzySmaller(5.0 - 0.0005, 5.0, 0.001));
}

TEST(CompareTest, isStrictFuzzySmallerBeyondEpsilon)
{
  EXPECT_TRUE(isStrictFuzzySmaller(5.0 - 0.00101, 5.0, 0.001));
  EXPECT_TRUE(isStrictFuzzySmaller(5.0, 6.0, 0.001));
  EXPECT_FALSE(isStrictFuzzySmaller(6.0, 5.0, 0.001));
  EXPECT_FALSE(isStrictFuzzySmaller(5.0 + 0.00101, 5.0, 0.001));
}
