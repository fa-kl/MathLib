/*****************************************************************************************
 * @file: conversion.cpp
 *
 * @brief: Conversion functions
 *
 * @details: This file provides various conversion functions.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include "conversion.hpp"

namespace mathlib
{

real_t rad2deg(const real_t rad)
{
  return rad * (180.0 / PI);
}

real_t deg2rad(const real_t deg)
{
  return deg * (PI / 180.0);
}

}  // namespace mathlib
