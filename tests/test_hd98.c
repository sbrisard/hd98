#include <glib.h>
#include <math.h>
#include <stdio.h>

#include "hd98/halm_dragon_1998.h"
#include "hd98/hd98.h"
#include "hd98/hooke.h"

void setup_hooke_tests();

void setup_halm_dragon_1998_tests();

void assert_array_equal(size_t size, double const *actual,
                        double const *expected, double rtol, double atol) {
  for (size_t i = 0; i < size; i++) {
    double const err = fabs(actual[i] - expected[i]);
    double const tol = atol + rtol * fabs(expected[i]);
    g_assert_cmpfloat(err, <=, tol);
  }
}

static HD98_Material *hooke_new_default() {
  double const kappa = 76700.;
  double const mu = 41600.;
  return hd98_hooke_new(kappa - 2 * mu / HD98_DIM, mu);
}

static HD98_Material *halm_dragon_1998_new_default() {
  double const kappa = 60700.;
  double const mu = 31300.;
  double const lambda = kappa - 2 * mu / HD98_DIM;
  double const alpha = 16000.;
  double const beta = 31000.;
  double const k0 = 0.11;
  double const k1 = 2.2;
  return hd98_halm_dragon_1998_new(lambda, mu, alpha, beta, k0, k1);
}

static void test_global_update() {
  HD98_Material const
      *mat[] = {halm_dragon_1998_new_default(), hooke_new_default()};

  size_t const n = 10;
  size_t phase[n];
  size_t m = 0; /* Total number of internal variables */
  for (size_t i = 0; i < n; i++) {
    phase[i] = i % 2;
    m += mat[phase[i]]->type->niv;
  }

  double delta_eps[n * HD98_SYM], eps1[n * HD98_SYM], sig2_act[n * HD98_SYM];
  double sig2_exp[n * HD98_SYM];
  for (size_t i = 0; i < n * HD98_SYM; i++) {
    delta_eps[i] = 0.;
    eps1[i] = 0.;
    sig2_act[i] = 0.;
    sig2_exp[i] = 0.;
  }

  double omega1[m], omega2_act[m], omega2_exp[m];
  double C2_act[n * HD98_SYM * HD98_SYM], C2_exp[n * HD98_SYM * HD98_SYM];

  for (size_t i = 0; i < n; i++) {
    eps1[HD98_SYM * i] = 1.e-4 * i;
  }
  for (size_t i = 0; i < m; i++) {
    omega1[i] = ((double) i) / ((double) (m - 1)) * 0.4;
  }

  hd98_global_update(n, phase, mat, delta_eps, eps1, omega1, sig2_act,
                     omega2_act, C2_act);
  double *omega1_i = omega1;
  double *omega2_i = omega2_exp;
  for (size_t i = 0; i < n; i++) {
    HD98_Material const *mat_i = mat[phase[i]];
    mat_i->type->update(mat_i, delta_eps + HD98_SYM * i, eps1 + HD98_SYM * i,
                        omega1_i, sig2_exp + HD98_SYM * i, omega2_i,
                        C2_exp + HD98_SYM * HD98_SYM * i);
    omega1_i += mat_i->type->niv;
    omega2_i += mat_i->type->niv;
  }

  assert_array_equal(n * HD98_SYM, sig2_act, sig2_exp, 1e-15, 1e-15);
  assert_array_equal(m, omega2_act, omega2_exp, 1e-15, 1e-15);
  assert_array_equal(n * HD98_SYM * HD98_SYM, C2_act, C2_exp, 1e-15, 1e-15);
}

static void test_solve_polarization_plus() {
  double const atol = 1e-15;
  double const rtol = 1e-15;

  HD98_Material const *mat = halm_dragon_1998_new_default();

  double const mu0 = 10000;
  double const nu0 = 0.3;
  double const lambda0 = 2 * mu0 * nu0 / (1 - 2 * nu0);

  double const eps1[] = {0., 0., 0., 0., 0., 1e-3};
  double const omega1[] = {0.3};

  double const O[] = {0., 0., 0., 0., 0., 0.};
  double sig1[HD98_SYM];
  double unused[1];
  /* TODO this should be a call to a function `current_stress`. */
  mat->type->update(mat, O, eps1, omega1, sig1, unused, NULL);

  double const delta_eps[] = {1e-4, -2e-4, 3e-4, -4e-4, 5e-4, -6e-4};
  double sig2[HD98_SYM];
  double omega2[1];
  mat->type->update(mat, delta_eps, eps1, omega1, sig2, omega2, NULL);
  double tr_delta_eps = 0.;
  for (size_t i = 0; i < HD98_DIM; i++) {
    tr_delta_eps += delta_eps[i];
  }
  double delta_tau[HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) {
    delta_tau[i] = sig2[i] - sig1[i] + 2 * mu0 * delta_eps[i];
    if (i < HD98_DIM) {
      delta_tau[i] += lambda0 * tr_delta_eps;
    }
  }

  double delta_eps_act[HD98_SYM];
  int status = hd98_solve_polarization_plus(mat, lambda0, mu0, delta_tau, eps1,
                                            omega1, delta_eps_act);
  g_assert_cmpint(status, ==, 0);
  assert_array_equal(HD98_SYM, delta_eps_act, delta_eps, rtol, atol);
}

void setup_hd98_tests() {
  g_test_add_func("/hd98/hd98_global_update", test_global_update);
  g_test_add_func("/hd98/hd98_solve_polarization_plus",
                  test_solve_polarization_plus);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);
  setup_hooke_tests();
  setup_halm_dragon_1998_tests();
  setup_hd98_tests();
  return g_test_run();
}
