#include <array>
#include <iostream>

#include <catch2/catch.hpp>

#include "hd98/hooke.hpp"
#include "test_hd98.hpp"

static void test_current_state(hd98::Hooke const& mat) {
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
    assert_array_equal(hd98::sym, sig_act.data(), sig_exp.data(), 1e-15, 1e-15);
    eps[i] = 0.;
  }
}

static void test_update(hd98::Hooke const& mat) {
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
    assert_array_equal(hd98::sym, sig2_act.data(), sig2_exp.data(), 1e-15,
                       1e-15);
    assert_array_equal(hd98::sym * hd98::sym, C2_act.data(), C2_exp.data(),
                       1e-15, 1e-15);
    delta_eps[i] = 0.;
  }
}

TEST_CASE("Hooke") {
  double const mu = 1.2;
  double const nu = 0.3;
  double const lambda = 2 * mu * nu / (1 - 2 * nu);
  hd98::Hooke const mat{lambda, mu};

  SECTION("current_state") { test_current_state(mat); }
  SECTION("Update") { test_update(mat); }
}