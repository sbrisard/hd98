#include <glib.h>
#include <math.h>
#include <stdio.h>

#include "hd98.h"

void assert_array_equal(size_t size, double *actual, double *expected,
                        double rtol, double atol) {
  for (size_t i = 0; i < size; i++) {
    double const err = fabs(actual[i] - expected[i]);
    double const tol = atol + rtol * fabs(expected[i]);
    g_assert_cmpfloat(err, <=, tol);
  }
}

void test_hd98_proportional_strain(double const *eps_dot) {
  double atol = 1e-15;
  double rtol = 1e-15;
  HalmDragon1998 *mat = halm_dragon_1998_new_default();

  double tr_eps_dot = eps_dot[0] + eps_dot[1] + eps_dot[2];
  double C_eps_dot_h = mat->lambda * tr_eps_dot;
  double H_eps_dot_h = 2 * mat->alpha * tr_eps_dot;
  double C_eps_dot[SYM];
  double H_eps_dot[SYM];
  double eps_dot_H_eps_dot = 0.;
  for (size_t i = 0; i < SYM; i++) {
    C_eps_dot[i] = 2. * mat->mu * eps_dot[i];
    H_eps_dot[i] = 4. * mat->beta * eps_dot[i];
    if (i < DIM) {
      C_eps_dot[i] += C_eps_dot_h;
      H_eps_dot[i] += H_eps_dot_h;
    }
    eps_dot_H_eps_dot += eps_dot[i] * H_eps_dot[i];
  }

  double t_damage = sqrt(2. * mat->k0_sqrt2 / eps_dot_H_eps_dot);
  double t_max = 2.0 * t_damage;
  size_t num_increments = 10;

  double delta_t = t_max / (double)num_increments;
  double delta_xi = delta_t / t_damage;
  double delta_eps_p[SYM], delta_eps_m[SYM];
  for (size_t i = 0; i < SYM; i++) {
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
        omega_exp[step] =
            (xi <= 1.) ? 0. : ((xi * xi - 1.) * mat->k0_sqrt2 / mat->k1_sqrt2);
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
  double sig[SYM], sig_exp[SYM];
  for (size_t step = 0; step < num_steps; step++) {
    t += sign[step] * delta_t;
    double *const delta_eps = (sign[step] == 1) ? delta_eps_p : delta_eps_m;
    double omega_new;
    mat->update(mat, delta_eps, eps, &omega, sig, &omega_new, NULL);
    for (size_t i = 0; i < SYM; i++) {
      eps[i] += delta_eps[i];
    }
    omega = omega_new;

    for (size_t i = 0; i < SYM; i++) {
      sig_exp[i] = t * (C_eps_dot[i] - omega_exp[step] * H_eps_dot[i]);
    }
    assert_array_equal(SYM, sig, sig_exp, rtol, atol);
  }

  mat->free(mat);
}

void test_hd98_setup_tests() {
  double *eps1 = g_new(double, SYM);
  double *eps2 = g_new(double, SYM);
  for (size_t i = 0; i < SYM; i++) {
    eps1[i] = i < DIM ? 1. : 0.;
    eps2[i] = 0.;
  }
  eps2[SYM - 1] = 1.;
  g_test_add_data_func_full("/HalmDragon1998/strain-driven/hydrostatic", eps1,
                            test_hd98_proportional_strain, g_free);
  g_test_add_data_func_full("/HalmDragon1998/strain-driven/deviatoric", eps2,
                            test_hd98_proportional_strain, g_free);
}

void test_global_update() {
  Material **mat = malloc(2 * sizeof(Material *));
  mat[0] = halm_dragon_1998_new_default();
  double kappa1 = 76700.;
  double mu1 = 41600.;
  mat[1] = hooke_new(kappa1 - 2 * mu1 / (double)DIM, mu1);

  size_t n = 10;
  uint8_t *phase = malloc(n * sizeof(size_t));
  for (size_t i = 0; i < n; i++)
    phase[i] = i % 2;

  double *delta_eps = calloc(n * SYM, sizeof(double));
  double *eps1 = calloc(n * SYM, sizeof(double));
  double *omega1 = calloc(n, sizeof(double));
  double *sig2 = calloc(n * SYM, sizeof(double));
  double *omega2 = calloc(n, sizeof(double));
  double *C2 = calloc(n * SYM * SYM, sizeof(double));

  for (size_t i = 0; i < n; i++) {
    omega1[i] = 0.;
    omega2[i] = 0.;
  }

  global_update(n, delta_eps, eps1, omega1, phase, mat, sig2, omega2, C2);
  double *C2_ijk = C2;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < SYM; j++) {
      for (size_t k = 0; k < SYM; k++) {
        printf("%g\t", *C2_ijk);
        ++C2_ijk;
      }
      printf("\n");
    }
  }

  free(C2);
  free(omega2);
  free(sig2);
  free(omega1);
  free(eps1);
  free(delta_eps);
  free(phase);
  free(mat);
}

int main(int argc, char **argv) {
  /* test_hd98_proportional_strain(); */
  /* test_global_update(); */
  g_test_init(&argc, &argv, NULL);
  test_hd98_setup_tests();
  return g_test_run();
}
