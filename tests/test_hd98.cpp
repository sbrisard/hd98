#include <array>
#include <cmath>
#include <concepts>
#include <cstdio>
#include <iterator>
#include <memory>
#include <numeric>

#include <catch2/catch.hpp>

#include "hd98/composite.hpp"
#include "hd98/halm_dragon_1998.hpp"
#include "hd98/hd98.hpp"
#include "hd98/hooke.hpp"

using Tensor2 = std::array<double, hd98::sym>;
using Tensor4 = std::array<double, hd98::sym * hd98::sym>;

template <typename T>
auto scale(T value) {
  return [value](T x) { return value * x; };
}

template <typename T>
auto increment(T value) {
  return [value](T x) { return x + value; };
}

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

static void test_hooke_current_state(hd98::Hooke const& mat) {
  Tensor2 eps{};
  Tensor2 sig_act{};
  Tensor2 sig_exp{};
  for (size_t i = 0; i < hd98::sym; i++) {
    eps[i] = 1.;
    mat.current_state(eps.data(), nullptr, sig_act.data());
    std::fill(sig_exp.begin(), sig_exp.end(), 0.0);
    sig_exp[i] = 2 * mat.mu * eps[i];
    if (i < hd98::dim) {
      for (size_t j = 0; j < hd98::dim; j++) {
        sig_exp[j] += mat.lambda;
      }
    }
    assert_approx_equal(sig_act.cbegin(), sig_act.cend(), sig_exp.cbegin(),
                        1e-15, 1e-15);
    eps[i] = 0.;
  }
}

static void test_hooke_update(hd98::Hooke const& mat) {
  Tensor2 eps1{};
  Tensor2 delta_eps{};
  Tensor2 sig2_act{};
  Tensor2 sig2_exp{};
  Tensor4 C2_act{};
  Tensor4 C2_exp{};
  for (size_t i = 0, ij = 0; i < hd98::sym; i++) {
    for (size_t j = 0; j < hd98::sym; j++, ij++) {
      C2_exp[ij] = (i < hd98::dim) && (j < hd98::dim) ? mat.lambda : 0.;
      if (i == j) {
        C2_exp[ij] += 2 * mat.mu;
      }
    }
  }
  for (size_t i = 0; i < hd98::sym; i++) {
    delta_eps[i] = 1.;
    mat.update(delta_eps.data(), eps1.data(), nullptr, sig2_act.data(), nullptr,
               C2_act.data());
    mat.current_state(delta_eps.data(), nullptr, sig2_exp.data());
    assert_approx_equal(sig2_act.cbegin(), sig2_act.cend(), sig2_exp.cbegin(),
                        1e-15, 1e-15);
    assert_approx_equal(C2_act.cbegin(), C2_act.cend(), C2_exp.cbegin(), 1e-15,
                        1e-15);
    delta_eps[i] = 0.;
  }
}
void test_halm_dragon_1998_current_state(hd98::HalmDragon1998 const& mat) {
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

void test_halm_dragon_1998_update_proportional_strain(
    hd98::HalmDragon1998 const& mat, Tensor2 const eps_dot) {
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

void test_composite_update() {
  hd98::Composite composite{hooke_new_default(),
                            halm_dragon_1998_new_default()};
  constexpr size_t n = 10;
  using ScalarField = std::array<double, n>;
  using Tensor2Field = std::array<double, n * hd98::sym>;
  using Tensor4Field = std::array<double, n * hd98::sym * hd98::sym>;

  std::array<int, n> phase;

  ScalarField omega1{}, omega2_act{}, omega2_exp{};
  Tensor2Field delta_eps{}, eps1{}, sig2_exp{}, sig2_act{};
  Tensor4Field C2_exp{}, C2_act{};

  for (size_t i = 0; i < n; i++) {
    phase[i] = i % 2;
    eps1[hd98::sym * i] = 1.e-4 * i;
    omega1[i] = ((double)i) / (n - 1.0) * 0.4;
  }

  composite.update(phase.size(), phase.data(), delta_eps.data(), eps1.data(),
                   omega1.data(), sig2_act.data(), omega2_act.data(),
                   C2_act.data());
  for (size_t i = 0; i < n; i++) {
    if (phase[i] == 0) {
      composite.hooke.update(
          delta_eps.data() + hd98::sym * i, eps1.data() + hd98::sym * i,
          omega1.data() + i, sig2_exp.data() + hd98::sym * i,
          omega2_exp.data() + i, C2_exp.data() + hd98::sym * hd98::sym * i);
    } else {
      composite.halm_dragon_1998.update(
          delta_eps.data() + hd98::sym * i, eps1.data() + hd98::sym * i,
          omega1.data() + i, sig2_exp.data() + hd98::sym * i,
          omega2_exp.data() + i, C2_exp.data() + hd98::sym * hd98::sym * i);
    }
  }

  assert_approx_equal(sig2_act.cbegin(), sig2_act.cend(), sig2_exp.cbegin(),
                      1e-15, 1e-15);
  assert_approx_equal(omega2_act.cbegin(), omega2_act.cend(),
                      omega2_exp.cbegin(), 1e-15, 1e-15);
  assert_approx_equal(C2_act.cbegin(), C2_act.cend(), C2_exp.cbegin(), 1e-15,
                      1e-15);
}

TEST_CASE("Hooke") {
  double const mu = 1.2;
  double const nu = 0.3;
  double const lambda = 2 * mu * nu / (1 - 2 * nu);
  hd98::Hooke const mat{lambda, mu};

  SECTION("Current state") { test_hooke_current_state(mat); }
  SECTION("Update") { test_hooke_update(mat); }
}

TEST_CASE("HalmDragon1998") {
  double const kappa = 60700.;
  double const mu = 31300.;
  double const lambda = kappa - 2 * mu / hd98::dim;
  double const alpha = 16000.;
  double const beta = 31000.;
  double const k0 = 0.11;
  double const k1 = 2.2;

  hd98::HalmDragon1998 mat{
      lambda, mu, alpha, beta, k0, k1, hd98::tangent_stiffness};

  SECTION("Current state") { test_halm_dragon_1998_current_state(mat); }

  SECTION("Update 1") {
    Tensor2 eps1{};
    std::fill(eps1.begin(), eps1.begin() + hd98::dim, 1.);
    test_halm_dragon_1998_update_proportional_strain(mat, eps1);
  }

  SECTION("Update 2") {
    Tensor2 eps2{};
    eps2[hd98::sym - 1] = 1.;
    test_halm_dragon_1998_update_proportional_strain(mat, eps2);
  }
}

TEST_CASE("Composite") {
  SECTION("Update") { test_composite_update(); }
}