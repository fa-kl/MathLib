/*****************************************************************************************
 * @file: IndexOutOfRangeError.hpp
 *
 * @brief: An index out of range error class.
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

class IndexOutOfRangeError : public MathLibError {

public:
  explicit IndexOutOfRangeError()
      : IndexOutOfRangeError("Index out of range.") {}

  explicit IndexOutOfRangeError(std::string msg)
      : MathLibError("IndexOutOfRangeError", std::move(msg)) {}
};

} // namespace mathlib