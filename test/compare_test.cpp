#include "compare.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace mathlib;

// test isFuzzyEqual
TEST(CompareTest,
     isFuzzyEqual_Input_ValOne5_ValTwo5point1_Epsilon0point001_Expect_False) {
  bool result = isFuzzyEqual(5, 5.1, 0.001);
  EXPECT_FALSE(result);
}

TEST(CompareTest,
     isFuzzyEqual_Input_ValOne5_ValTwo5point001_Epsilon0point001_Expect_False) {
  bool result = isFuzzyEqual(5, 5.001, 0.001);
  EXPECT_FALSE(result);
}

TEST(CompareTest,
     isFuzzyEqual_Input_ValOne3_ValTwo3_Epsilon0point001_Expect_True) {
  bool result = isFuzzyEqual(3, 3, 0.001);
  EXPECT_TRUE(result);
}

TEST(
    CompareTest,
    isFuzzyEqual_Input_ValOne5_ValTwo5point00099_Epsilon0point001_Expect_True) {
  bool result = isFuzzyEqual(5, 5.00099, 0.001);
  EXPECT_TRUE(result);
}

// test isFuzzyGreater
TEST(
    CompareTest,
    isFuzzyGreater_Input_ValOne5_ValTwo5point001_Epsilon0point001_Expect_False) {
  bool result = isFuzzyGreater(5, 5.001, 0.001);
  EXPECT_FALSE(result);
}

TEST(CompareTest,
     isFuzzyGreater_Input_ValOne5_ValTwo5_Epsilon0point001_Expect_True) {
  bool result = isFuzzyEqual(5, 5, 0.001);
  EXPECT_TRUE(result);
}

// test isFuzzySmaller
TEST(CompareTest,
     isFuzzySmaller_Input_ValOne5_ValTwo6_Epsilon0point001_Expect_True) {
  bool result = isFuzzySmaller(5, 6, 0.001);
  EXPECT_TRUE(result);
}

TEST(CompareTest,
     isFuzzySmaller_Input_ValOne5_ValTwo4_Epsilon0point001_Expect_False) {
  bool result = isFuzzySmaller(5, 4, 0.001);
  EXPECT_FALSE(result);
}

// test isStrictFuzzyGreater
TEST(CompareTest,
     isStrictFuzzyGreater_Input_ValOne6_ValTwo5_Epsilon0point001_Expect_True) {
  bool result = isStrictFuzzyGreater(6, 5, 0.001);
  EXPECT_TRUE(result);
}

TEST(CompareTest,
     isStrictFuzzyGreater_Input_ValOne5_ValTwo6_Epsilon0point001_Expect_False) {
  bool result = isStrictFuzzyGreater(5, 6, 0.001);
  EXPECT_FALSE(result);
}

// test isStrictFuzzySmaller
TEST(CompareTest,
     isStrictFuzzyGreater_Input_ValOne6_ValTwo5_Epsilon0point001_Expect_False) {
  bool result = isStrictFuzzySmaller(6, 5, 0.001);
  EXPECT_FALSE(result);
}

TEST(CompareTest,
     isStrictFuzzyGreater_Input_ValOne5_ValTwo6_Epsilon0point001_Expect_True) {
  bool result = isStrictFuzzySmaller(5, 6, 0.001);
  EXPECT_TRUE(result);
}
