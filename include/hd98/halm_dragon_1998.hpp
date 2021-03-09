#pragma once

#include <numbers>
#include <sstream>
#include "hd98/hd98.hpp"

namespace hd98 {
class HalmDragon1998 {
 public:
  double const lambda;
  double const mu;
  double const
      alpha; /* Should be > 0 (opposite convention to the paper of Xianda. */
  double const beta;
  double const k0_sqrt2;
  double const k1_sqrt2;
  int const stiffness_type;

  HalmDragon1998(double lambda, double mu, double alpha, double beta, double k0,
                 double k1, int stiffness_type)
      : lambda(lambda),
        mu(mu),
        alpha(alpha),
        beta(beta),
        k0_sqrt2(k0 * std::numbers::sqrt2),
        k1_sqrt2(k1 * std::numbers::sqrt2),
        stiffness_type(stiffness_type) {}

  [[nodiscard]] std::string repr() const {
    std::ostringstream stream;
    stream << "Hooke"
           << "{lambda=" << lambda << ",mu=" << mu << ",alpha=" << alpha
           << ",beta=" << beta << ",k0=" << 0.5 * k0_sqrt2 * std::numbers::sqrt2
           << ",k1=" << 0.5 * k1_sqrt2 * std::numbers::sqrt2 << "}";
    return stream.str();
  }

  void current_state(double const *eps, double const *omega,
                     double *sig) const {
    double two_mu = 2 * mu - 4 * beta * omega[0];
    double lambda_tr_eps = 0.;
    for (size_t i = 0; i < dim; i++) {
      lambda_tr_eps += eps[i];
    }
    lambda_tr_eps *= lambda - 2 * alpha * omega[0];
    for (size_t i = 0; i < dim; i++) {
      sig[i] = lambda_tr_eps + two_mu * eps[i];
    }
    for (size_t i = dim; i < sym; i++) {
      sig[i] = two_mu * eps[i];
    }
  }

  void update(double const *delta_eps, double const *eps1, double const *omega1,
              double *sig2, double *omega2, double *C2) const {
    double eps2[sym];
    for (size_t i = 0; i < sym; i++) eps2[i] = eps1[i] + delta_eps[i];

    double tr_eps2 = 0.;
    for (size_t i = 0; i < dim; i++) tr_eps2 += eps2[i];

    double two_alpha_tr_eps2 = 2. * alpha * tr_eps2;
    double four_beta = 4. * beta;
    double H_eps2[sym];
    for (size_t i = 0; i < sym; i++) H_eps2[i] = four_beta * eps2[i];
    for (size_t i = 0; i < dim; i++) H_eps2[i] += two_alpha_tr_eps2;

    double eps2_H_eps2 = 0.;
    for (size_t i = 0; i < sym; i++) eps2_H_eps2 += eps2[i] * H_eps2[i];

    double f_tr = 0.5 * eps2_H_eps2 - (k0_sqrt2 + k1_sqrt2 * omega1[0]);

    double delta_omega = f_tr > 0. ? f_tr / k1_sqrt2 : 0.;
    omega2[0] = omega1[0] + delta_omega;

    double lambda_tr_eps2 = lambda * tr_eps2;
    double two_mu = 2. * mu;
    for (size_t i = 0; i < sym; i++)
      sig2[i] = two_mu * eps2[i] - omega2[0] * H_eps2[i];
    for (size_t i = 0; i < dim; i++) sig2[i] += lambda_tr_eps2;

    if (C2 != nullptr) {
      double lambda_sec = lambda - 2. * omega2[0] * alpha;
      double two_mu_sec = 2. * (mu - 2. * omega2[0] * beta);
      double *C2_ij = C2;
      if ((stiffness_type == tangent_stiffness) && (f_tr > 0)) {
        for (size_t i = 0; i < sym; i++) {
          double aux = -H_eps2[i] / k1_sqrt2;
          for (size_t j = 0; j < sym; j++) {
            *C2_ij = aux * H_eps2[j];
            if (i == j) *C2_ij += two_mu_sec;
            if ((i < dim) && (j < dim)) *C2_ij += lambda_sec;
            ++C2_ij;
          }
        }
      } else {
        for (size_t i = 0; i < sym; i++) {
          for (size_t j = 0; j < sym; j++) {
            *C2_ij = 0;
            if (i == j) *C2_ij += two_mu_sec;
            if ((i < dim) && (j < dim)) *C2_ij += lambda_sec;
            ++C2_ij;
          }
        }
      }
    }
  }
};

std::ostream &operator<<(std::ostream &os, const HalmDragon1998 &mat) {
  return os << mat.repr();
}

}  // namespace hd98