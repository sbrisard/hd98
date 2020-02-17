#include <stdlib.h>
#include <string.h>

#include "hd98/hooke.h"

void hd98_hooke_free(HD98_Material *mat) {
  HD98_HookeData *data = mat->data;
  free(data->C);
  free(data);
  free(mat);
}

void hd98_hooke_update(HD98_Material const *mat, double const *delta_eps,
                       double const *eps1, double const *unused1, double *sig2,
                       double *unused2, double *C2) {
  HD98_HookeData *data = mat->data;
  double eps2[HD98_SYM];
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

HD98_MaterialType const HD98_Hooke = {.name = "Hooke",
                                      .num_int_var = 0,
                                      .free = hd98_hooke_free,
                                      .update = hd98_hooke_update};

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
