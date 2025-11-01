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

#include <exception>
#include <string>
#include <utility>

#include "MathLibError.hpp"
#include "config.hpp"

namespace mathlib
{

class IndexOutOfRangeError : public MathLibError
{
public:
  explicit IndexOutOfRangeError() : IndexOutOfRangeError("Index out of range.") {}

  explicit IndexOutOfRangeError(const size_t i) : IndexOutOfRangeError("Index " + std::to_string(i) + " out of range.") {}

  explicit IndexOutOfRangeError(const index_t i) : IndexOutOfRangeError("Index " + std::to_string(i) + " out of range.") {}

  explicit IndexOutOfRangeError(std::string msg) : MathLibError("IndexOutOfRangeError", std::move(msg)) {}
};

}  // namespace mathlib