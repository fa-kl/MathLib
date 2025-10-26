/*****************************************************************************************
 * @file: MathLibError.hpp
 *
 * @brief: A general MathLib error class.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#pragma once

#include <exception>
#include <string>
#include <utility>

namespace mathlib {

class MathLibError : public std::exception {
protected:
  std::string m_type;
  std::string m_msg;

public:
  explicit MathLibError(std::string msg)
      : m_type("MathLibError"), m_msg(std::move(msg)) {}

  explicit MathLibError(std::string type, std::string msg)
      : m_type(std::move(type)), m_msg(std::move(msg)) {}

  [[nodiscard]] const char *what() const noexcept override {
    static std::string str = "[" + m_type + "]: " + m_msg;
    return str.c_str();
  }
};

} // namespace mathlib