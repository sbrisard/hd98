#include <glib.h>

#include "hd98/hooke.h"
#include "test_hd98.h"

static void test_material_type() {
  g_assert_cmpstr(HD98_Hooke.name, ==, "Hooke");
  g_assert_cmpuint(HD98_Hooke.niv, ==, 0);
}

static void test_new() {
  double mu = 1.2;
  double nu = 0.3;
  double lambda = 2 * mu * nu / (1 - 2 * nu);
  HD98_Material *mat = hd98_hooke_new(lambda, mu);
  g_assert_true(mat->type == &HD98_Hooke);
  HD98_HookeData *data = mat->data;
  g_assert_cmpfloat(data->lambda, ==, lambda);
  g_assert_cmpfloat(data->mu, ==, mu);
}

static HD98_Material *hooke_new_default() {
  double mu = 1.2;
  double nu = 0.3;
  double lambda = 2 * mu * nu / (1 - 2 * nu);
  return hd98_hooke_new(lambda, mu);
}

static void test_current_state() {
  HD98_Material *mat = hooke_new_default();
  HD98_HookeData *data = mat->data;
  double eps[HD98_SYM], sig_act[HD98_SYM], sig_exp[HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) {
    eps[i] = 0.;
  }
  for (size_t i = 0; i < HD98_SYM; i++) {
    eps[i] = 1.;
    mat->type->current_state(mat, eps, NULL, sig_act);
    for (size_t j = 0; j < HD98_SYM; j++) {
      sig_exp[j] = 0.;
    }
    sig_exp[i] = 2 * data->mu * eps[i];
    if (i < HD98_DIM) {
      for (size_t j = 0; j < HD98_DIM; j++) {
        sig_exp[j] += data->lambda;
      }
    }
    assert_array_equal(HD98_SYM, sig_act, sig_exp, 1e-15, 1e-15);
    eps[i] = 0.;
  }
}

static void test_update() {
  HD98_Material *mat = hooke_new_default();
  HD98_HookeData *data = mat->data;
  double eps1[HD98_SYM];
  double delta_eps[HD98_SYM];
  double sig2_act[HD98_SYM];
  double sig2_exp[HD98_SYM];
  double C2_act[HD98_SYM * HD98_SYM];
  double C2_exp[HD98_SYM * HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) {
    eps1[i] = 0.;
    delta_eps[i] = 0.;
  }
  for (size_t i = 0, ij = 0; i < HD98_SYM; i++) {
    for (size_t j = 0; j < HD98_SYM; j++, ij++) {
      C2_exp[ij] = (i < HD98_DIM) && (j < HD98_DIM) ? data->lambda : 0.;
      if (i == j) {
        C2_exp[ij] += 2 * data->mu;
      }
    }
  }
  for (size_t i = 0; i < HD98_SYM; i++) {
    delta_eps[i] = 1.;
    mat->type->update(mat, delta_eps, eps1, NULL, sig2_act, NULL, C2_act);
    mat->type->current_state(mat, delta_eps, NULL, sig2_exp);
    assert_array_equal(HD98_SYM, sig2_act, sig2_exp, 1e-15, 1e-15);
    assert_array_equal(HD98_SYM * HD98_SYM, C2_act, C2_exp, 1e-15, 1e-15);
    delta_eps[i] = 0.;
  }
}

void setup_hooke_tests() {
  g_test_add_func("/Hooke/HD98_MaterialType", test_material_type);
  g_test_add_func("/Hooke/new", test_new);
  g_test_add_func("/Hooke/current_state", test_current_state);
  g_test_add_func("/Hooke/update", test_update);
}
