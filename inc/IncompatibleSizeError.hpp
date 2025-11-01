/*****************************************************************************************
 * @file: IncompatibleSizeError.hpp
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

namespace mathlib
{

class IncompatibleSizeError : public MathLibError
{
public:
  explicit IncompatibleSizeError() : IncompatibleSizeError("Sizes are not compatible for this operatation.") {}

  IncompatibleSizeError(const size_t rows1, const size_t cols1, const size_t rows2, const size_t cols2)
      : IncompatibleSizeError()
  {
    m_msg += " (";
    m_msg += std::to_string(rows1);
    m_msg += ", ";
    m_msg += std::to_string(cols1);
    m_msg += ") != (";
    m_msg += std::to_string(rows2);
    m_msg += ", ";
    m_msg += std::to_string(cols2);
    m_msg += ")";
  }

  explicit IncompatibleSizeError(std::string msg) : MathLibError("IncompatibleSizeError", std::move(msg)) {}
};

}  // namespace mathlib