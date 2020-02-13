#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "hd98/hd98.h"

typedef struct HD98_HookeData_ {
  double lambda;
  double mu;
  double *C;
} HD98_HookeData;

void hd98_hooke_free(HD98_Material *mat) {
  HD98_HookeData *data = mat->data;
  free(data->C);
  free(data);
  free(mat);
}

void hd98_hooke_update(HD98_Material const *mat, double const *delta_eps,
                       double const *eps1, double const *unused1, double *sig2,
                       double *unused2, double *C2) {
  double eps2[HD98_SYM];
  HD98_HookeData *data = mat->data;
  for (size_t i = 0; i < HD98_SYM; i++) eps2[i] = eps1[i] + delta_eps[i];

  double lambda_tr_eps2 = 0.;
  for (size_t i = 0; i < HD98_DIM; i++) lambda_tr_eps2 += eps2[i];
  lambda_tr_eps2 *= data->lambda;

  double two_mu = 2. * data->mu;
  for (size_t i = 0; i < HD98_SYM; i++) sig2[i] = two_mu * eps2[i];
  for (size_t i = 0; i < HD98_DIM; i++) sig2[i] += lambda_tr_eps2;
  if (C2) {
    memcpy(C2, data->C, HD98_SYM * HD98_SYM * sizeof(double));
  }
}

HD98_MaterialType const HD98_Hooke = {
    .name = "Hooke", .free = hd98_hooke_free, .update = hd98_hooke_update};

HD98_Material *hd98_hooke_new(double lambda, double mu) {
  HD98_HookeData *data = malloc(sizeof(HD98_HookeData));
  data->lambda = lambda;
  data->mu = mu;
  data->C = malloc(HD98_SYM * HD98_SYM * sizeof(double));
  double *C_ij = data->C;
  for (int i = 0; i < HD98_SYM; i++) {
    for (int j = 0; j < HD98_SYM; j++) {
      *C_ij = 0.;
      if (i == j) *C_ij += 2. * mu;
      if ((i < HD98_DIM) && (j < HD98_DIM)) *C_ij += lambda;
      ++C_ij;
    }
  }
  HD98_Material *hooke = malloc(sizeof(HD98_Material));
  hooke->type = &HD98_Hooke;
  hooke->data = data;
  return hooke;
}

typedef struct HD98_HalmDragon1998Data_ {
  double lambda;
  double mu;
  double alpha; /* Should be > 0 (opposite convention to Xianda's paper. */
  double beta;
  double k0_sqrt2;
  double k1_sqrt2;
} HD98_HalmDragon1998Data;

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

HD98_Material *hd98_halm_dragon_1998_new_default() {
  double kappa = 60700.;
  double mu = 31300.;
  double alpha = 16000.;
  double beta = 31000.;
  double k0 = 0.11;
  double k1 = 2.2;
  return hd98_halm_dragon_1998_new(kappa - 2 * mu / HD98_DIM, mu, alpha, beta,
                                   k0, k1);
}

void hd98_global_update(size_t n, double const *delta_eps, double const *eps1,
                        double const *omega1, size_t const *phase,
                        HD98_Material **mat, double *sig2, double *omega2,
                        double *C2) {
  double const *delta_eps_i = delta_eps;
  double const *eps1_i = eps1;
  double const *omega1_i = omega1;
  HD98_Material const *mat_i;
  double *sig2_i = sig2;
  double *omega2_i = omega2;
  double *C2_i = C2;

  for (size_t i = 0; i < n; i++) {
    mat_i = mat[phase[i]];
    mat_i->type->update(mat_i, delta_eps_i, eps1_i, omega1_i, sig2_i, omega2_i,
                        C2_i);

    delta_eps_i += HD98_SYM;
    eps1_i += HD98_SYM;
    ++omega1_i; /* This assumes that there is only one internal variable. */
    sig2_i += HD98_SYM;
    ++omega2_i;
    C2_i += HD98_SYM * HD98_SYM;
  }
}
