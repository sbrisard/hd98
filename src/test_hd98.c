#include <math.h>
#include <glib.h>

#include "hd98.h"

void test_hd98_proportional_strain(double const *eps_star) {
  double const atol = 1e-15;
  double const rtol = 1e-15;
  HalmDragon1998 *mat = halm_dragon_1998_new_default();

  double const tr_eps_star = eps_star[0] + eps_star[1] + eps_star[2];
  double const C_eps_star_h = mat->lambda * tr_eps_star;
  double const H_eps_star_h = 2 * mat->alpha * tr_eps_star;
  double C_eps_star[SYM];
  double H_eps_star[SYM];
  double eps_star_H_eps_star = 0.;
  for (size_t i = 0; i < SYM; i++) {
    C_eps_star[i] = 2. * mat->mu * eps_star[i];
    H_eps_star[i] = 4. * mat->beta * eps_star[i];
    if (i < DIM) {
      C_eps_star[i] += C_eps_star_h;
      H_eps_star[i] += H_eps_star_h;
    }
    eps_star_H_eps_star += eps_star[i] * H_eps_star[i];
  }

  double const xi_damage = sqrt(2. * mat->k0_sqrt2 / eps_star_H_eps_star);
  double const xi_max = 2.0 * xi_damage;
  size_t const num_steps = 11;
  double const delta_xi = xi_max / (num_steps - 1.0);
  double delta_eps[SYM];
  for (size_t i = 0; i < SYM; i++) {
    delta_eps[i] = delta_xi * eps_star[i];
  }

  double omega = 0.;
  double eps[] = {0., 0., 0., 0., 0., 0.};
  double sig[SYM];
  for (size_t k = 1; k < num_steps; k++) {
    double omega_new;
    mat->update(mat, delta_eps, eps, &omega, sig, &omega_new, NULL);
    for (size_t i = 0; i < SYM; i++) {
      eps[i] += delta_eps[i];
    }
    omega = omega_new;

    double const xi = k * delta_xi;
    double const xi_red = xi / xi_damage;
    double const omega_exp =
        (xi_red <= 1.)
            ? 0.
            : ((xi_red * xi_red - 1.) * mat->k0_sqrt2 / mat->k1_sqrt2);
    for (size_t i = 0; i < SYM; i++) {
      double const sig_exp_i = xi * (C_eps_star[i] - omega_exp * H_eps_star[i]);
      double const res = fabs(sig[i] - sig_exp_i);
      double const tol = rtol * fabs(sig_exp_i) + atol;
      g_assert_cmpfloat(res, <=, tol);
    }
  }

  /* Unloading */

  mat->free(mat);
}

void test_hd98_setup_tests() {
  double *eps1 = g_new(double, SYM);
  for (size_t i = 0; i < SYM; i++) {
    eps1[i] = i < DIM ? 1. : 0.;
  }
  g_test_add_data_func_full("/HalmDragon1998/strain-driven/hydrostatic", eps1,
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
