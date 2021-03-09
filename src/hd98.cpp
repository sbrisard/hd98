#include "hd98/hd98.hpp"

namespace hd98 {

// int hd98_solve_polarization_plus(HD98_Material const *mat, double lambda0,
//                                 double mu0, double const *delta_tau,
//                                 double const *eps1, double const *iv1,
//                                 double *delta_eps) {
//  /* TODO These values should not be hard-coded. */
//  double atol = 1e-15;
//  double rtol = 1e-15;
//  size_t max_iter = 10;
//
//  double iv2[mat->type->niv], sig1[sym], sig2[sym],
//      C2[sym * sym];
//
//  /* A: matrix of NR iterations; b: residual; x: correction to delta_eps */
//  gsl_matrix *A = gsl_matrix_calloc(sym, sym);
//  gsl_vector *x = gsl_vector_calloc(sym);
//  gsl_vector *b = gsl_vector_calloc(sym);
//  gsl_permutation *p = gsl_permutation_alloc(sym);
//
//  /* Compute initial stress */
//  /* TODO this should be a call to a function `current_stress`. */
//  for (size_t i = 0; i < sym; i++) delta_eps[i] = 0.;
//  mat->type->update(mat, delta_eps, eps1, iv1, sig1, iv2, NULL);
//
//  /* Define iter outside the loop in order to be able to return an
//     error code. */
//  size_t iter = 0;
//  for (; iter <= max_iter; iter++) {
//    mat->type->update(mat, delta_eps, eps1, iv1, sig2, iv2, C2);
//    /* Compute matrix */
//    for (size_t i = 0, ij = 0; i < sym; i++) {
//      for (size_t j = 0; j < sym; j++, ij++) {
//        double A_ij = C2[ij];
//        if (i == j) A_ij += 2 * mu0;
//        if ((i < dim) && (j < dim)) A_ij += lambda0;
//        gsl_matrix_set(A, i, j, A_ij);
//      }
//    }
//    /* Compute residual */
//    double tr_delta_eps = 0;
//    for (size_t i = 0; i < dim; i++) tr_delta_eps += delta_eps[i];
//    bool converged = true;
//    for (size_t i = 0; i < sym; i++) {
//      double b_i = delta_tau[i] - (sig2[i] - sig1[i]) - 2 * mu0 *
//      delta_eps[i]; if (i < dim) b_i -= lambda0 * tr_delta_eps; converged
//      = converged && (fabs(b_i) <= rtol * fabs(delta_tau[i]) + atol);
//      gsl_vector_set(b, i, b_i);
//    }
//    if (converged) break;
//    /* Compute correction */
//    int s;
//    gsl_linalg_LU_decomp(A, p, &s);
//    gsl_linalg_LU_solve(A, p, b, x);
//    /* Apply correction */
//    for (size_t i = 0; i < sym; i++) {
//      delta_eps[i] += gsl_vector_get(x, i);
//    }
//  }
//
//  gsl_matrix_free(A);
//  gsl_vector_free(b);
//  gsl_vector_free(x);
//  gsl_permutation_free(p);
//
//  return iter > max_iter;
//}
//
// int hd98_solve_polarizations_plus(size_t n, size_t const *phase,
//                                  HD98_Material const **mat, double lambda0,
//                                  double mu0, double const *delta_tau,
//                                  double const *eps1, double const *iv1,
//                                  double *delta_eps) {
//  double const *delta_tau_i = delta_tau;
//  double const *eps1_i = eps1;
//  double const *iv1_i = iv1;
//  HD98_Material const *mat_i;
//  double *delta_eps_i = delta_eps;
//
//  for (size_t i = 0; i < n; i++) {
//    mat_i = mat[phase[i]];
//    int err = hd98_solve_polarization_plus(mat_i, lambda0, mu0, delta_tau_i,
//                                           eps1_i, iv1_i, delta_eps_i);
//    if (err) {
//        return i;
//    }
//    delta_tau_i += sym;
//    eps1_i += sym;
//    iv1_i += mat_i->type->niv;
//    delta_eps_i += sym;
//  }
//  return 0;
//}
}  // namespace hd98