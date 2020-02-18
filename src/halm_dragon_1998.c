#include <math.h>
#include <stdlib.h>

#include "hd98/halm_dragon_1998.h"

void hd98_halm_dragon_1998_free(HD98_Material *mat) {
  free(mat->data);
  free(mat);
}

void hd98_halm_dragon_1998_update(HD98_Material const *mat,
                                  double const *delta_eps, double const *eps1,
                                  double const *omega1, double *sig2,
                                  double *omega2, double *C2) {
  HD98_HalmDragon1998Data *data = mat->data;
  double eps2[HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) eps2[i] = eps1[i] + delta_eps[i];

  double tr_eps2 = 0.;
  for (size_t i = 0; i < HD98_DIM; i++) tr_eps2 += eps2[i];

  double two_alpha_tr_eps2 = 2. * data->alpha * tr_eps2;
  double four_beta = 4. * data->beta;
  double H_eps2[HD98_SYM];
  for (size_t i = 0; i < HD98_SYM; i++) H_eps2[i] = four_beta * eps2[i];
  for (size_t i = 0; i < HD98_DIM; i++) H_eps2[i] += two_alpha_tr_eps2;

  double eps2_H_eps2 = 0.;
  for (size_t i = 0; i < HD98_SYM; i++) eps2_H_eps2 += eps2[i] * H_eps2[i];

  double f_tr =
      0.5 * eps2_H_eps2 - (data->k0_sqrt2 + data->k1_sqrt2 * omega1[0]);

  double delta_omega = f_tr > 0. ? f_tr / data->k1_sqrt2 : 0.;
  omega2[0] = omega1[0] + delta_omega;

  double lambda_tr_eps2 = data->lambda * tr_eps2;
  double two_mu = 2. * data->mu;
  for (size_t i = 0; i < HD98_SYM; i++)
    sig2[i] = two_mu * eps2[i] - omega2[0] * H_eps2[i];
  for (size_t i = 0; i < HD98_DIM; i++) sig2[i] += lambda_tr_eps2;

  if (C2 != NULL) {
    double lambda_sec = data->lambda - 2. * omega2[0] * data->alpha;
    double two_mu_sec = 2. * (data->mu - 2. * omega2[0] * data->beta);
    double *C2_ij = C2;
    if (f_tr > 0) {
      for (size_t i = 0; i < HD98_SYM; i++) {
        double aux = -H_eps2[i] / data->k1_sqrt2;
        for (size_t j = 0; j < HD98_SYM; j++) {
          *C2_ij = aux * H_eps2[j];
          if (i == j) *C2_ij += two_mu_sec;
          if ((i < HD98_DIM) && (j < HD98_DIM)) *C2_ij += lambda_sec;
          ++C2_ij;
        }
      }
    } else {
      for (size_t i = 0; i < HD98_SYM; i++) {
        for (size_t j = 0; j < HD98_SYM; j++) {
          *C2_ij = 0;
          if (i == j) *C2_ij += two_mu_sec;
          if ((i < HD98_DIM) && (j < HD98_DIM)) *C2_ij += lambda_sec;
          ++C2_ij;
        }
      }
    }
  }
}

HD98_MaterialType const HD98_HalmDragon1998 = {
    .name = "HalmDragon1998",
    .num_int_var = 1,
    .free = hd98_halm_dragon_1998_free,
    .update = hd98_halm_dragon_1998_update};

HD98_Material *hd98_halm_dragon_1998_new(double lambda, double mu, double alpha,
                                         double beta, double k0, double k1) {
  HD98_HalmDragon1998Data *data = malloc(sizeof(HD98_HalmDragon1998Data));
  data->lambda = lambda;
  data->mu = mu;
  data->alpha = alpha;
  data->beta = beta;
  data->k0_sqrt2 = k0 * M_SQRT2;
  data->k1_sqrt2 = k1 * M_SQRT2;
  HD98_Material *mat = malloc(sizeof(HD98_Material));
  mat->type = &HD98_HalmDragon1998;
  mat->data = data;
  return mat;
}