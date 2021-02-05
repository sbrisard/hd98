#include <cstdlib>
#include <cstring>

#include "hd98/hooke.hpp"

static void hooke_free(HD98_Material *mat) {
  auto data = static_cast<HD98_HookeData *>(mat->data);
  free(data->C);
  free(data);
  free(mat);
}

static void hooke_current_state(HD98_Material const *mat, double const *eps,
                                double const *unused, double *sig) {
  auto data = static_cast<HD98_HookeData *>(mat->data);
  double lambda_tr_eps = 0.;
  for (size_t i = 0; i < HD98_DIM; i++) {
    lambda_tr_eps += eps[i];
  }
  lambda_tr_eps *= data->lambda;
  double two_mu = 2 * data->mu;
  for (size_t i = 0; i < HD98_DIM; i++) {
    sig[i] = lambda_tr_eps + two_mu * eps[i];
  }
  for (size_t i = HD98_DIM; i < HD98_SYM; i++) {
    sig[i] = two_mu * eps[i];
  }
}

static void hooke_update(HD98_Material const *mat, double const *delta_eps,
                         double const *eps1, double const *unused1,
                         double *sig2, double *unused2, double *C2) {
  auto data = static_cast<HD98_HookeData *>(mat->data);
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

HD98_MaterialType const HD98_Hooke{"Hooke", 0, hooke_free, hooke_current_state,
                                   hooke_update};

HD98_Material *hd98_hooke_new(double lambda, double mu) {
  auto data = static_cast<HD98_HookeData *>(malloc(sizeof(HD98_HookeData)));
  data->lambda = lambda;
  data->mu = mu;
  data->C = static_cast<double *>(malloc(HD98_SYM * HD98_SYM * sizeof(double)));
  double *C_ij = data->C;
  for (int i = 0; i < HD98_SYM; i++) {
    for (int j = 0; j < HD98_SYM; j++) {
      *C_ij = 0.;
      if (i == j) *C_ij += 2. * mu;
      if ((i < HD98_DIM) && (j < HD98_DIM)) *C_ij += lambda;
      ++C_ij;
    }
  }
  auto hooke = static_cast<HD98_Material *>(malloc(sizeof(HD98_Material)));
  hooke->type = &HD98_Hooke;
  hooke->data = data;
  return hooke;
}
