#include "conversion.hpp"

#include <iostream>

#include <gtest/gtest.h>

using namespace mathlib;

TEST(ConversionTest, Rad2Deg)
{
  EXPECT_NEAR(rad2deg(M_PI / 2), 90.0, 1e-6);
  EXPECT_NEAR(rad2deg(M_PI), 180.0, 1e-6);
  EXPECT_NEAR(rad2deg(3 * M_PI / 2), 270.0, 1e-6);
  EXPECT_NEAR(rad2deg(2 * M_PI), 360.0, 1e-6);
  EXPECT_NEAR(rad2deg(0.0), 0.0, 1e-6);
  EXPECT_NEAR(rad2deg(-M_PI / 2), -90.0, 1e-6);
  EXPECT_NEAR(rad2deg(0.001), 0.057295779513082320876798154814105, 1e-6);
}

TEST(ConversionTest, Deg2Rad)
{
  EXPECT_NEAR(deg2rad(90.0), M_PI / 2, 1e-6);
  EXPECT_NEAR(deg2rad(180.0), M_PI, 1e-6);
  EXPECT_NEAR(deg2rad(270.0), 3 * M_PI / 2, 1e-6);
  EXPECT_NEAR(deg2rad(360.0), 2 * M_PI, 1e-6);
  EXPECT_NEAR(deg2rad(0.0), 0.0, 1e-6);
  EXPECT_NEAR(deg2rad(-90.0), -M_PI / 2, 1e-6);
  EXPECT_NEAR(deg2rad(0.057295779513082320876798154814105), 0.001, 1e-6);
}