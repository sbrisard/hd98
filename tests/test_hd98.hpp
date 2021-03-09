#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <iterator>

#include "hd98/halm_dragon_1998.hpp"
#include "hd98/hd98.hpp"
#include "hd98/hooke.hpp"
#include "hd98/composite.hpp"

using Tensor2 = std::array<double, hd98::sym>;
using Tensor4 = std::array<double, hd98::sym * hd98::sym>;

inline void assert_approx_equal(double act, double exp, double rtol,
                                double atol) {
  INFO("exp = " << exp << ", act = " << act);
  REQUIRE(abs(act - exp) <= rtol * abs(exp) + atol);
}

template <typename It>
requires std::input_iterator<It>&&
    std::same_as<std::iter_value_t<It>, double> void
    assert_approx_equal(It act_start, It act_end, It exp_start, double rtol,
                        double atol) {
  for (auto act = act_start, exp = exp_start; act != act_end; ++act, ++exp) {
    assert_approx_equal(*act, *exp, rtol, atol);
  }
}

hd98::Hooke hooke_new_default() {
  double kappa = 76700.;
  double mu = 41600.;
  return hd98::Hooke{kappa - 2 * mu / hd98::dim, mu};
}

hd98::HalmDragon1998 halm_dragon_1998_new_default() {
  double kappa = 60700.;
  double mu = 31300.;
  double lambda = kappa - 2 * mu / hd98::dim;
  double alpha = 16000.;
  double beta = 31000.;
  double k0 = 0.11;
  double k1 = 2.2;
  return hd98::HalmDragon1998{
      lambda, mu, alpha, beta, k0, k1, hd98::tangent_stiffness};
}
