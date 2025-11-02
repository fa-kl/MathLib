/*****************************************************************************************
 * @file: main.cpp
 *
 * @brief: This is the main file.
 *
 * @details: This file may be used to try some functionality.
 *
 * @author: fakl
 * @date: October 2025
 *
 ****************************************************************************************/

#include <iostream>

#include "mathlib.hpp"

using namespace mathlib;

/*****************************************************************************************
 *  Main Function
 ****************************************************************************************/
int main(void)
{
  /* TODO: try it out yourself! */

  Vector v = {1.0, 2.0, 3.0};
  std::cout << "v:\n" << v << "\n";
  std::cout << "mean(v) = " << mean(v) << "\n";
  std::cout << "std(v) = " << mathlib::std(v) << "\n";
  std::cout << "var(v) = " << var(v) << "\n";

  std::cout << "\n\n";

  Matrix A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  std::cout << "A:\n" << A << "\n\n";
  std::cout << "mean(A) =\n" << mean(A, 2) << "\n\n";
  std::cout << "std(A) =\n" << mathlib::std(A, 2) << "\n\n";
  std::cout << "var(A) =\n" << var(A, 2) << "\n";

  return 0;
}
