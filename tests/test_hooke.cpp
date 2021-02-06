#include <cstdio>
#include "hd98/hooke.hpp"
#include "test_hd98.hpp"

static void test_current_state(hd98::Hooke const&mat) {
  printf("Hooke/test_current_state...");
  double eps[hd98::sym], sig_act[hd98::sym], sig_exp[hd98::sym];
  for (size_t i = 0; i < hd98::sym; i++) {
    eps[i] = 0.;
  }
  for (size_t i = 0; i < hd98::sym; i++) {
    eps[i] = 1.;
    mat.current_state(eps, NULL, sig_act);
    for (size_t j = 0; j < hd98::sym; j++) {
      sig_exp[j] = 0.;
    }
    sig_exp[i] = 2 * mat.mu * eps[i];
    if (i < hd98::dim) {
      for (size_t j = 0; j < hd98::dim; j++) {
        sig_exp[j] += mat.lambda;
      }
    }
    assert_array_equal(hd98::sym, sig_act, sig_exp, 1e-15, 1e-15);
    eps[i] = 0.;
  }
  printf(" OK\n");
}

static void test_update(hd98::Hooke const& mat) {
  printf("Hooke/test_update...");
  double eps1[hd98::sym];
  double delta_eps[hd98::sym];
  double sig2_act[hd98::sym];
  double sig2_exp[hd98::sym];
  double C2_act[hd98::sym * hd98::sym];
  double C2_exp[hd98::sym * hd98::sym];
  for (size_t i = 0; i < hd98::sym; i++) {
    eps1[i] = 0.;
    delta_eps[i] = 0.;
  }
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
    mat.update(delta_eps, eps1, NULL, sig2_act, NULL, C2_act);
    mat.current_state(delta_eps, NULL, sig2_exp);
    assert_array_equal(hd98::sym, sig2_act, sig2_exp, 1e-15, 1e-15);
    assert_array_equal(hd98::sym * hd98::sym, C2_act, C2_exp, 1e-15, 1e-15);
    delta_eps[i] = 0.;
  }
  printf(" OK\n");
}

void setup_hooke_tests() {
  double const mu = 1.2;
  double const nu = 0.3;
  double const lambda = 2 * mu * nu / (1 - 2 * nu);
  hd98::Hooke const mat{lambda, mu};
  test_current_state(mat);
  test_update(mat);
}
