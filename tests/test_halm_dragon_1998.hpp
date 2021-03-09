#pragma once

#include <numeric>
#include <vector>

#include <catch2/catch.hpp>

#include "hd98/halm_dragon_1998.hpp"
#include "hd98/hooke.hpp"
#include "test_hd98.hpp"

namespace test_halm_dragon_1998 {
template <typename T>
auto scale(T value) {
  return [value](T x) { return value * x; };
}

template <typename T>
auto increment(T value) {
  return [value](T x) { return x + value; };
}

static void test_current_state(hd98::HalmDragon1998 const &mat) {
  Tensor2 eps{1.2, -3.4, 5.6, -7.8, 9., -10.11};
  std::array<double, 1> omega{0.4};
  Tensor2 sig_act{};
  Tensor2 sig_exp{};
  mat.current_state(eps.data(), omega.data(), sig_act.data());

  hd98::Hooke hooke{mat.lambda - 2 * omega[0] * mat.alpha,
                    mat.mu - 2 * omega[0] * mat.beta};
  hooke.current_state(eps.data(), nullptr, sig_exp.data());
  assert_approx_equal(sig_act.cbegin(), sig_act.cend(), sig_exp.cbegin(), 1e-15,
                      1e-15);
}

static void test_update_proportional_strain(hd98::HalmDragon1998 const &mat,
                                            Tensor2 const eps_dot) {
  double atol = 1e-15;
  double rtol = 1e-15;

  auto tr_eps_dot =
      std::accumulate(eps_dot.cbegin(), eps_dot.cbegin() + hd98::dim, 0.0);

  Tensor2 C_eps_dot, H_eps_dot;
  std::transform(eps_dot.cbegin(), eps_dot.cend(), C_eps_dot.begin(),
                 scale(2 * mat.mu));
  std::transform(eps_dot.cbegin(), eps_dot.cend(), H_eps_dot.begin(),
                 scale(4 * mat.beta));
  std::transform(C_eps_dot.cbegin(), C_eps_dot.cbegin() + hd98::dim,
                 C_eps_dot.begin(), increment(mat.lambda * tr_eps_dot));
  std::transform(H_eps_dot.cbegin(), H_eps_dot.cbegin() + hd98::dim,
                 H_eps_dot.begin(), increment(2 * mat.alpha * tr_eps_dot));
  auto eps_dot_H_eps_dot = std::inner_product(eps_dot.cbegin(), eps_dot.cend(),
                                              H_eps_dot.cbegin(), 0.0);

  double t_damage = sqrt(2 * mat.k0_sqrt2 / eps_dot_H_eps_dot);
  double t_max = 2 * t_damage;
  size_t num_increments = 10;

  double delta_t = t_max / (double)num_increments;
  double delta_xi = delta_t / t_damage;
  Tensor2 delta_eps_p, delta_eps_m;
  std::transform(eps_dot.cbegin(), eps_dot.cend(), delta_eps_p.begin(),
                 scale(delta_t));
  std::transform(eps_dot.cbegin(), eps_dot.cend(), delta_eps_m.begin(),
                 scale(-delta_t));

  size_t num_steps = num_increments * (num_increments + 1);
  std::vector<int> sign(num_steps);
  std::vector<double> omega_exp(num_steps);
  for (size_t increment = 1, step = 0; increment <= num_increments;
       increment++) {
    for (size_t i = 1; i <= increment; i++, step++) {
      sign[step] = 1;
      if (i == increment) {
        double xi = increment * delta_xi;
        omega_exp[step] =
            (xi <= 1.) ? 0. : ((xi * xi - 1.) * mat.k0_sqrt2 / mat.k1_sqrt2);
      } else {
        omega_exp[step] = omega_exp[step - 1];
      }
    }
    for (size_t i = 1; i <= increment; i++, step++) {
      sign[step] = -1;
      omega_exp[step] = omega_exp[step - 1];
    }
  }

  double t = 0.;
  double omega = 0.;
  Tensor2 eps{}, sig, sig_exp;
  for (size_t step = 0; step < num_steps; step++) {
    t += sign[step] * delta_t;
    auto delta_eps = (sign[step] == 1) ? delta_eps_p : delta_eps_m;
    double omega_new;
    mat.update(delta_eps.data(), eps.data(), &omega, sig.data(), &omega_new,
               nullptr);
    std::transform(eps.cbegin(), eps.cend(), delta_eps.cbegin(), eps.begin(),
                   std::plus());
    omega = omega_new;

    for (size_t i = 0; i < hd98::sym; i++) {
      sig_exp[i] = t * (C_eps_dot[i] - omega_exp[step] * H_eps_dot[i]);
    }
    assert_approx_equal(sig.cbegin(), sig.cend(), sig_exp.cbegin(), rtol, atol);
  }
}
}  // namespace test_halm_dragon_1998

