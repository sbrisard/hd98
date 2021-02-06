#include <iostream>
#include <numeric>
#include <vector>

#include "hd98/halm_dragon_1998.hpp"
#include "hd98/hooke.hpp"
#include "test_hd98.hpp"

static void test_current_state(hd98::HalmDragon1998 const &mat) {
  std::cout << "HalmDragon1998/test_current_state...";
  Tensor2 eps{1.2, -3.4, 5.6, -7.8, 9., -10.11};
  std::array<double, 1> omega{0.4};
  Tensor2 sig_act{};
  Tensor2 sig_exp{};
  mat.current_state(eps.data(), omega.data(), sig_act.data());

  hd98::Hooke hooke{mat.lambda - 2 * omega[0] * mat.alpha,
                    mat.mu - 2 * omega[0] * mat.beta};
  hooke.current_state(eps.data(), nullptr, sig_exp.data());
  assert_array_equal(hd98::sym, sig_act.data(), sig_exp.data(), 1e-15, 1e-15);
  std::cout << " OK\n";
}

static void test_update_proportional_strain(hd98::HalmDragon1998 const &mat,
                                            Tensor2 const eps_dot) {
  std::cout << "HalmDragon1998/test_update_proportional_strain...";
  double atol = 1e-15;
  double rtol = 1e-15;

  auto tr_eps_dot =
      std::accumulate(eps_dot.cbegin(), eps_dot.cbegin() + hd98::dim, 0.0);

  double C_eps_dot_h = mat.lambda * tr_eps_dot;
  double H_eps_dot_h = 2 * mat.alpha * tr_eps_dot;
  double C_eps_dot[hd98::sym];
  double H_eps_dot[hd98::sym];
  double eps_dot_H_eps_dot = 0.;
  for (size_t i = 0; i < hd98::sym; i++) {
    C_eps_dot[i] = 2 * mat.mu * eps_dot[i];
    H_eps_dot[i] = 4 * mat.beta * eps_dot[i];
    if (i < hd98::dim) {
      C_eps_dot[i] += C_eps_dot_h;
      H_eps_dot[i] += H_eps_dot_h;
    }
    eps_dot_H_eps_dot += eps_dot[i] * H_eps_dot[i];
  }

  double t_damage = sqrt(2 * mat.k0_sqrt2 / eps_dot_H_eps_dot);
  double t_max = 2 * t_damage;
  size_t num_increments = 10;

  double delta_t = t_max / (double)num_increments;
  double delta_xi = delta_t / t_damage;
  double delta_eps_p[hd98::sym], delta_eps_m[hd98::sym];
  for (size_t i = 0; i < hd98::sym; i++) {
    delta_eps_p[i] = delta_t * eps_dot[i];
    delta_eps_m[i] = -delta_eps_p[i];
  }

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
  double eps[] = {0., 0., 0., 0., 0., 0.};
  double sig[hd98::sym], sig_exp[hd98::sym];
  for (size_t step = 0; step < num_steps; step++) {
    t += sign[step] * delta_t;
    double *delta_eps = (sign[step] == 1) ? delta_eps_p : delta_eps_m;
    double omega_new;
    mat.update(delta_eps, eps, &omega, sig, &omega_new, nullptr);
    for (size_t i = 0; i < hd98::sym; i++) {
      eps[i] += delta_eps[i];
    }
    omega = omega_new;

    for (size_t i = 0; i < hd98::sym; i++) {
      sig_exp[i] = t * (C_eps_dot[i] - omega_exp[step] * H_eps_dot[i]);
    }
    assert_array_equal(hd98::sym, sig, sig_exp, rtol, atol);
  }
  std::cout << " OK\n";
}

void setup_halm_dragon_1998_tests() {
  double const kappa = 60700.;
  double const mu = 31300.;
  double const lambda = kappa - 2 * mu / hd98::dim;
  double const alpha = 16000.;
  double const beta = 31000.;
  double const k0 = 0.11;
  double const k1 = 2.2;

  hd98::HalmDragon1998 mat{
      lambda, mu, alpha, beta, k0, k1, hd98::tangent_stiffness};

  Tensor2 eps1{};
  Tensor2 eps2{};
  std::fill(eps1.begin(), eps1.begin()+hd98::dim, 1.);
  eps2[hd98::sym - 1] = 1.;
  test_current_state(mat);
  test_update_proportional_strain(mat, eps1);
  test_update_proportional_strain(mat, eps2);
}
