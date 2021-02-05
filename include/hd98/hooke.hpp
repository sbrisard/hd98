#pragma once

#include <array>
#include <cstring>

#include "hd98/hd98.hpp"

namespace hd98 {
DllExport std::array<double, HD98_SYM * HD98_SYM> stiffness_matrix(
    double lambda, double mu);

class Hooke {
 public:
  double const lambda;
  double const mu;
  std::array<double, HD98_SYM * HD98_SYM> const C;

  Hooke(double lambda, double mu)
      : lambda(lambda), mu(mu), C(stiffness_matrix(lambda, mu)) {}

  void current_state(double const *eps, double const *unused,
                     double *sig) const {
    double lambda_tr_eps = 0.;
    for (size_t i = 0; i < HD98_DIM; i++) {
      lambda_tr_eps += eps[i];
    }
    lambda_tr_eps *= lambda;
    double two_mu = 2 * mu;
    for (size_t i = 0; i < HD98_DIM; i++) {
      sig[i] = lambda_tr_eps + two_mu * eps[i];
    }
    for (size_t i = HD98_DIM; i < HD98_SYM; i++) {
      sig[i] = two_mu * eps[i];
    }
  }

  void update(double const *delta_eps, double const *eps1,
              double const *unused1, double *sig2, double *unused2,
              double *C2) const {
    double eps2[HD98_SYM];
    for (size_t i = 0; i < HD98_SYM; i++) eps2[i] = eps1[i] + delta_eps[i];

    double lambda_tr_eps2 = 0.;
    for (size_t i = 0; i < HD98_DIM; i++) lambda_tr_eps2 += eps2[i];
    lambda_tr_eps2 *= lambda;

    double two_mu = 2. * mu;
    for (size_t i = 0; i < HD98_SYM; i++) sig2[i] = two_mu * eps2[i];
    for (size_t i = 0; i < HD98_DIM; i++) sig2[i] += lambda_tr_eps2;
    if (C2) {
      std::memcpy(C2, C.data(), HD98_SYM * HD98_SYM * sizeof(double));
    }
  }

 private:
  static std::array<double, HD98_SYM * HD98_SYM> stiffness_matrix(double lambda,
                                                                  double mu) {
    std::array<double, HD98_SYM * HD98_SYM> C{};
    for (size_t i = 0, ij = 0; i < HD98_SYM; i++) {
      for (size_t j = 0; j < HD98_SYM; j++, ij++) {
        C[ij] = 0.0;
        if (i == j) C[ij] += 2. * mu;
        if ((i < HD98_DIM) && (j < HD98_DIM)) C[ij] += lambda;
      }
    }
    return C;
  }
};
}  // namespace hd98