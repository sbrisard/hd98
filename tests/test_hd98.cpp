#include <cstdio>
#include <memory>

#include <catch2/catch.hpp>

#include "test_hd98.hpp"
#include "test_hooke.hpp"
#include "test_halm_dragon_1998.hpp"

static void test_global_update() {
  hd98::Composite composite{hooke_new_default(),
                            halm_dragon_1998_new_default()};
  constexpr size_t n = 10;
  using ScalarField = std::array<double, n>;
  using Tensor2Field = std::array<double, n * hd98::sym>;
  using Tensor4Field = std::array<double, n * hd98::sym * hd98::sym>;

  std::array<size_t, n> phase;

  ScalarField omega1{}, omega2_act{}, omega2_exp{};
  Tensor2Field delta_eps{}, eps1{}, sig2_exp{}, sig2_act{};
  Tensor4Field C2_exp{}, C2_act{};

  for (size_t i = 0; i < n; i++) {
    phase[i] = i % 2;
    eps1[hd98::sym * i] = 1.e-4 * i;
    omega1[i] = ((double)i) / (n - 1.0) * 0.4;
  }

  composite.update(phase, delta_eps.data(), eps1.data(), omega1.data(),
                   sig2_act.data(), omega2_act.data(), C2_act.data());
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

// static void test_solve_polarization_plus() {
//  double atol = 1e-15;
//  double rtol = 1e-15;
//
//  HD98_Material *mat = halm_dragon_1998_new_default();
//
//  double mu0 = 10000;
//  double nu0 = 0.3;
//  double lambda0 = 2 * mu0 * nu0 / (1 - 2 * nu0);
//
//  double eps1[] = {0., 0., 0., 0., 0., 1e-3};
//  double omega1[] = {0.3};
//
//  double O[] = {0., 0., 0., 0., 0., 0.};
//  double sig1[sym];
//  double unused[1];
//  /* TODO this should be a call to a function `current_stress`. */
//  mat->type->update(mat, O, eps1, omega1, sig1, unused, NULL);
//
//  double delta_eps[] = {1e-4, -2e-4, 3e-4, -4e-4, 5e-4, -6e-4};
//  double sig2[sym];
//  double omega2[1];
//  mat->type->update(mat, delta_eps, eps1, omega1, sig2, omega2, NULL);
//  double tr_delta_eps = 0.;
//  for (size_t i = 0; i < dim; i++) {
//    tr_delta_eps += delta_eps[i];
//  }
//  double delta_tau[sym];
//  for (size_t i = 0; i < sym; i++) {
//    delta_tau[i] = sig2[i] - sig1[i] + 2 * mu0 * delta_eps[i];
//    if (i < dim) {
//      delta_tau[i] += lambda0 * tr_delta_eps;
//    }
//  }
//
//  double delta_eps_act[sym];
//  int status = hd98_solve_polarization_plus(mat, lambda0, mu0, delta_tau,
//  eps1,
//                                            omega1, delta_eps_act);
//  g_assert_cmpint(status, ==, 0);
//  assert_array_equal(sym, delta_eps_act, delta_eps, rtol, atol);
//}

TEST_CASE("Hooke") {
  double const mu = 1.2;
  double const nu = 0.3;
  double const lambda = 2 * mu * nu / (1 - 2 * nu);
  hd98::Hooke const mat{lambda, mu};

  SECTION("current_state") { test_hooke::test_current_state(mat); }
  SECTION("Update") { test_hooke::test_update(mat); }
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

  SECTION("current_state") { test_halm_dragon_1998::test_current_state(mat); }

  SECTION("update 1") {
    Tensor2 eps1{};
    std::fill(eps1.begin(), eps1.begin() + hd98::dim, 1.);
    test_halm_dragon_1998::test_update_proportional_strain(mat, eps1);
  }

  SECTION("update 2") {
    Tensor2 eps2{};
    eps2[hd98::sym - 1] = 1.;
    test_halm_dragon_1998::test_update_proportional_strain(mat, eps2);
  }
}

TEST_CASE("hd98") {
  SECTION("global_update") { test_global_update(); }
}