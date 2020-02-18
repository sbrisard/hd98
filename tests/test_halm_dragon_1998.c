#include <glib.h>
#include <math.h>

#include "hd98/halm_dragon_1998.h"
#include "test_hd98.h"

static void test_material_type() {
  g_assert_cmpstr(HD98_HalmDragon1998.name, ==, "HalmDragon1998");
  g_assert_cmpuint(HD98_HalmDragon1998.num_int_var, ==, 1);
}

static void test_new() {
  double const kappa = 60700.;
  double const mu = 31300.;
  double const lambda = kappa - 2 * mu / HD98_DIM;
  double const alpha = 16000.;
  double const beta = 31000.;
  double const k0 = 0.11;
  double const k1 = 2.2;
  HD98_Material const *mat =
      hd98_halm_dragon_1998_new(lambda, mu, alpha, beta, k0, k1);
  g_assert_true(mat->type == &HD98_HalmDragon1998);
  HD98_HalmDragon1998Data *data = mat->data;
  g_assert_cmpfloat(data->lambda, ==, lambda);
  g_assert_cmpfloat(data->mu, ==, mu);
  g_assert_cmpfloat(data->alpha, ==, alpha);
  g_assert_cmpfloat(data->beta, ==, beta);
  g_assert_cmpfloat(data->k0_sqrt2, ==, k0 * M_SQRT2);
  g_assert_cmpfloat(data->k1_sqrt2, ==, k1 * M_SQRT2);
}

static void test_update_proportional_strain(gconstpointer data) {
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

  HD98_Material *mat =
      hd98_halm_dragon_1998_new(lambda, mu, alpha, beta, k0, k1);

  double tr_eps_dot = 0.;
  for (size_t i = 0; i < HD98_DIM; i++) {
    tr_eps_dot += eps_dot[i];
  }
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

void setup_halm_dragon_1998_tests() {
  double *eps1 = g_new(double, HD98_SYM);
  double *eps2 = g_new(double, HD98_SYM);
  for (size_t i = 0; i < HD98_SYM; i++) {
    eps1[i] = i < HD98_DIM ? 1. : 0.;
    eps2[i] = 0.;
  }
  eps2[HD98_SYM - 1] = 1.;
  g_test_add_func("/HalmDragon1998/HD98_MaterialType", test_material_type);
  g_test_add_func("/HalmDragon1998/new", test_new);
  g_test_add_data_func_full("/HalmDragon1998/update/hydrostatic", eps1,
                            test_update_proportional_strain, g_free);
  g_test_add_data_func_full("/HalmDragon1998/update/deviatoric", eps2,
                            test_update_proportional_strain, g_free);
}