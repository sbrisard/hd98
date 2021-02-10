#pragma once

#include <array>
#include <cstring>
#include <sstream>

#include "hd98/hd98.hpp"

namespace hd98 {
class Hooke : Material {
 public:
  double const lambda;
  double const mu;
  std::array<double, sym * sym> const C;

  Hooke(double lambda, double mu)
      : lambda(lambda), mu(mu), C(stiffness_matrix(lambda, mu)) {}

  [[nodiscard]] std::string repr() const override {
    std::ostringstream stream;
    stream << "Hooke"
           << "{lambda=" << lambda << ","
           << "mu=" << mu << "}";
    return stream.str();
  }

  void current_state(double const *eps, double const *unused,
                     double *sig) const override {
    double lambda_tr_eps = 0.;
    for (size_t i = 0; i < dim; i++) {
      lambda_tr_eps += eps[i];
    }
    lambda_tr_eps *= lambda;
    double two_mu = 2 * mu;
    for (size_t i = 0; i < dim; i++) {
      sig[i] = lambda_tr_eps + two_mu * eps[i];
    }
    for (size_t i = dim; i < sym; i++) {
      sig[i] = two_mu * eps[i];
    }
  }

  void update(double const *delta_eps, double const *eps1,
              double const *unused1, double *sig2, double *unused2,
              double *C2) const override {
    double eps2[sym];
    for (size_t i = 0; i < sym; i++) eps2[i] = eps1[i] + delta_eps[i];

    double lambda_tr_eps2 = 0.;
    for (size_t i = 0; i < dim; i++) lambda_tr_eps2 += eps2[i];
    lambda_tr_eps2 *= lambda;

    double two_mu = 2. * mu;
    for (size_t i = 0; i < sym; i++) sig2[i] = two_mu * eps2[i];
    for (size_t i = 0; i < dim; i++) sig2[i] += lambda_tr_eps2;
    if (C2) {
      std::memcpy(C2, C.data(), sym * sym * sizeof(double));
    }
  }

 private:
  static std::array<double, sym * sym> stiffness_matrix(double lambda,
                                                        double mu) {
    std::array<double, sym * sym> C{};
    for (size_t i = 0, ij = 0; i < sym; i++) {
      for (size_t j = 0; j < sym; j++, ij++) {
        C[ij] = 0.0;
        if (i == j) C[ij] += 2. * mu;
        if ((i < dim) && (j < dim)) C[ij] += lambda;
      }
    }
    return C;
  }
};
}  // namespace hd98