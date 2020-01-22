#include <glib.h>
#include <math.h>
#include <stdio.h>

#include "hd98/hd98.h"

void assert_array_equal(size_t size, double *actual, double *expected,
                        double rtol, double atol) {
  for (size_t i = 0; i < size; i++) {
    double const err = fabs(actual[i] - expected[i]);
    double const tol = atol + rtol * fabs(expected[i]);
    g_assert_cmpfloat(err, <=, tol);
  }
}

void test_hd98_proportional_strain(gconstpointer data) {
  double const *eps_dot = data;
  double atol = 1e-15;
  double rtol = 1e-15;

  double kappa = 60700.;
  double mu = 31300.;
  double lambda = kappa - 2 * mu / HD98_DIM;
  double alpha = 16000.;
  double beta = 31000.;
  double k0 = 0.11;
  double k1 = 2.2;

  Material *mat = hd98_halm_dragon_1998_new_default();

  double tr_eps_dot = eps_dot[0] + eps_dot[1] + eps_dot[2];
  double C_eps_dot_h = lambda * tr_eps_dot;
  double H_eps_dot_h = 2 * alpha * tr_eps_dot;
  double C_eps_dot[HD98_SYM];
  double H_eps_dot[HD98_SYM];
  double eps_dot_H_eps_dot = 0.;
  for (size_t i = 0; i < HD98_SYM; i++) {
    C_eps_dot[i] = 2. * mu * eps_dot[i];
    H_eps_dot[i] = 4. * beta * eps_dot[i];
    if (i < HD98_DIM) {
      C_eps_dot[i] += C_eps_dot_h;
      H_eps_dot[i] += H_eps_dot_h;
    }
    eps_dot_H_eps_dot += eps_dot[i] * H_eps_dot[i];
  }

  double t_damage = sqrt(2. * k0 * M_SQRT2 / eps_dot_H_eps_dot);
  double t_max = 2.0 * t_damage;
  size_t num_increments = 10;

  double delta_t = t_max / (double)num_increments;
  double delta_xi = delta_t / t_damage;
  double delta_eps_p[HD98_SYM], delta_eps_m[HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) {
    delta_eps_p[i] = delta_t * eps_dot[i];
    delta_eps_m[i] = -delta_eps_p[i];
  }

  size_t num_steps = num_increments * (num_increments + 1);
  int sign[num_steps];
  double omega_exp[num_steps];
  for (size_t increment = 1, step = 0; increment <= num_increments;
       increment++) {
    for (size_t i = 1; i <= increment; i++, step++) {
      sign[step] = 1;
      if (i == increment) {
        double const xi = increment * delta_xi;
        omega_exp[step] = (xi <= 1.) ? 0. : ((xi * xi - 1.) * k0 / k1);
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
  double sig[HD98_SYM], sig_exp[HD98_SYM];
  for (size_t step = 0; step < num_steps; step++) {
    t += sign[step] * delta_t;
    double *const delta_eps = (sign[step] == 1) ? delta_eps_p : delta_eps_m;
    double omega_new;
    mat->type->update(mat, delta_eps, eps, &omega, sig, &omega_new, NULL);
    for (size_t i = 0; i < HD98_SYM; i++) {
      eps[i] += delta_eps[i];
    }
    omega = omega_new;

    for (size_t i = 0; i < HD98_SYM; i++) {
      sig_exp[i] = t * (C_eps_dot[i] - omega_exp[step] * H_eps_dot[i]);
    }
    assert_array_equal(HD98_SYM, sig, sig_exp, rtol, atol);
  }

  mat->type->free(mat);
}

void test_hd98_global_update() {
  Material **mat = malloc(2 * sizeof(Material *));
  mat[0] = hd98_halm_dragon_1998_new_default();
  double kappa1 = 76700.;
  double mu1 = 41600.;
  mat[1] = hd98_hooke_new(kappa1 - 2 * mu1 / (double) HD98_DIM, mu1);

  size_t n = 10;
  uint8_t *phase = malloc(n * sizeof(size_t));
  for (size_t i = 0; i < n; i++) phase[i] = (uint8_t)(i % 2);

  double *delta_eps = calloc(n * HD98_SYM, sizeof(double));
  double *eps1 = calloc(n * HD98_SYM, sizeof(double));
  double *omega1 = calloc(n, sizeof(double));

  double *sig2_act = calloc(n * HD98_SYM, sizeof(double));
  double *omega2_act = calloc(n, sizeof(double));
  double *C2_act = calloc(n * HD98_SYM * HD98_SYM, sizeof(double));

  double *sig2_exp = calloc(n * HD98_SYM, sizeof(double));
  double *omega2_exp = calloc(n, sizeof(double));
  double *C2_exp = calloc(n * HD98_SYM * HD98_SYM, sizeof(double));

  for (size_t i = 0; i < n; i++) {
    eps1[HD98_SYM * i] = 1.e-4 * i;
    omega1[i] = ((double)i) / ((double)(n - 1)) * 0.4;
  }

  hd98_global_update(n,
                     delta_eps,
                     eps1,
                     omega1,
                     phase,
                     mat,
                     sig2_act,
                     omega2_act,
                     C2_act);
  for (size_t i = 0; i < n; i++) {
    Material *mat_i = mat[phase[i]];
    mat_i->type->update(mat_i, delta_eps + HD98_SYM * i, eps1 + HD98_SYM * i,
                        omega1 + i, sig2_exp + HD98_SYM * i, omega2_exp + i,
                        C2_exp + HD98_SYM * HD98_SYM * i);
  }

  assert_array_equal(n * HD98_SYM, sig2_act, sig2_exp, 1e-15, 1e-15);
  assert_array_equal(n, omega2_act, omega2_exp, 1e-15, 1e-15);
  assert_array_equal(n * HD98_SYM * HD98_SYM, C2_act, C2_exp, 1e-15, 1e-15);

  free(C2_exp);
  free(omega2_exp);
  free(sig2_exp);
  free(C2_act);
  free(omega2_act);
  free(sig2_act);
  free(omega1);
  free(eps1);
  free(delta_eps);
  free(phase);
  free(mat);
}

void test_hd98_setup_tests() {
  double *eps1 = g_new(double, HD98_SYM);
  double *eps2 = g_new(double, HD98_SYM);
  for (size_t i = 0; i < HD98_SYM; i++) {
    eps1[i] = i < HD98_DIM ? 1. : 0.;
    eps2[i] = 0.;
  }
  eps2[HD98_SYM - 1] = 1.;
  g_test_add_data_func_full("/HalmDragon1998/strain-driven/hydrostatic", eps1,
                            test_hd98_proportional_strain, g_free);
  g_test_add_data_func_full("/HalmDragon1998/strain-driven/deviatoric", eps2,
                            test_hd98_proportional_strain, g_free);
  g_test_add_func("/HalmDragon1998/hd98_global_update", test_hd98_global_update);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);
  test_hd98_setup_tests();
  return g_test_run();
}
