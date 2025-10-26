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

#include "MathLibError.hpp"
#include <exception>
#include <string>
#include <utility>

namespace mathlib {

class DivByZeroError : public MathLibError {

public:
  explicit DivByZeroError() : DivByZeroError("A division by zero occured.") {}

  explicit DivByZeroError(std::string msg)
      : MathLibError("DivByZeroError", std::move(msg)) {}
};

} // namespace mathlib