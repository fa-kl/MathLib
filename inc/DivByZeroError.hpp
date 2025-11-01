/*****************************************************************************************
 * @file: DivByZeroError.hpp
 *
 * @brief: A division by zero error class.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <exception>
#include <string>
#include <utility>

#include "MathLibError.hpp"

namespace mathlib
{

class DivByZeroError : public MathLibError
{
public:
  explicit DivByZeroError() : DivByZeroError("A division by zero occured.") {}

  explicit DivByZeroError(std::string msg) : MathLibError("DivByZeroError", std::move(msg)) {}
};

}  // namespace mathlib