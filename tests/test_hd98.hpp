#pragma once

#include <array>
#include <cmath>

#include "hd98/hd98.hpp"

using Tensor2 = std::array<double, hd98::sym>;
using Tensor4 = std::array<double, hd98::sym * hd98::sym>;

inline void assert_approx_equal(double act, double exp, double rtol,
                                double atol) {
  INFO("exp = " << exp << ", act = " << act);
  REQUIRE(abs(act - exp) <= rtol * abs(exp) + atol);
}

template <typename It>
void assert_approx_equal(It act_start, It act_end, It exp_start, double rtol,
                         double atol) {
  for (auto act = act_start, exp = exp_start; act != act_end; ++act, ++exp) {
    assert_approx_equal(*act, *exp, rtol, atol);
  }
}