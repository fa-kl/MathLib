#include "conversion.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace mathlib;
TEST(ConversionTest, Rad2Deg) {
  // Test basic conversion
  EXPECT_NEAR(rad2deg(M_PI / 2), 90.0, 1e-6);      // 90 degrees
  EXPECT_NEAR(rad2deg(M_PI), 180.0, 1e-6);         // 180 degrees
  EXPECT_NEAR(rad2deg(3 * M_PI / 2), 270.0, 1e-6); // 270 degrees
  EXPECT_NEAR(rad2deg(2 * M_PI), 360.0, 1e-6);     // 360 degrees

  // Test zero and negative angles
  EXPECT_NEAR(rad2deg(0.0), 0.0, 1e-6);         // Zero degrees
  EXPECT_NEAR(rad2deg(-M_PI / 2), -90.0, 1e-6); // Negative 90 degrees

  // Test small angles
  EXPECT_NEAR(rad2deg(0.001), 0.057295779513082320876798154814105, 1e-6);
}

TEST(ConversionTest, Deg2Rad) {
  // Test basic conversion
  EXPECT_NEAR(deg2rad(90.0), M_PI / 2, 1e-6);      // π/2 radians
  EXPECT_NEAR(deg2rad(180.0), M_PI, 1e-6);         // π radians
  EXPECT_NEAR(deg2rad(270.0), 3 * M_PI / 2, 1e-6); // 3π/2 radians
  EXPECT_NEAR(deg2rad(360.0), 2 * M_PI, 1e-6);     // 2π radians

  // Test zero and negative angles
  EXPECT_NEAR(deg2rad(0.0), 0.0, 1e-6);         // Zero radians
  EXPECT_NEAR(deg2rad(-90.0), -M_PI / 2, 1e-6); // Negative π/2 radians

  // Test small angles
  EXPECT_NEAR(deg2rad(0.057295779513082320876798154814105), 0.001, 1e-6);
}