#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIM 3
#define SYM 6

typedef struct Material Material;

typedef void material_free_t(Material*);
typedef void material_update_t(Material*, double*, double*, double*, double*,
                               double*, double*);

struct Material {
  material_free_t* free;
  material_update_t* update;
};

typedef struct Hooke {
  struct Material;
  double lambda;
  double mu;
  double* C;
} Hooke;

void hooke_free(Hooke* mat) {
  free(mat->C);
  free(mat);
}

void hooke_update(Hooke* mat, double* delta_eps, double* eps1, double* unused1,
                  double* sig2, double* unused2, double* C2) {
  double eps2[SYM];
  for (size_t i = 0; i < SYM; i++) eps2[i] = eps1[i] + delta_eps[i];

  double lambda_tr_eps2 = 0.;
  for (size_t i = 0; i < DIM; i++) lambda_tr_eps2 += eps2[i];
  lambda_tr_eps2 *= mat->lambda;

  double two_mu = 2. * mat->mu;
  for (size_t i = 0; i < SYM; i++) sig2[i] = two_mu * eps2[i];
  for (size_t i = 0; i < DIM; i++) sig2[i] += lambda_tr_eps2;
  if (C2) {
    memcpy(C2, mat->C, SYM * SYM * sizeof(double));
  }
}

Hooke* hooke_new(double lambda, double mu) {
  Hooke* hooke = malloc(sizeof(Hooke));
  hooke->lambda = lambda;
  hooke->mu = mu;
  hooke->C = malloc(SYM * SYM * sizeof(double));
  double* C_ij = hooke->C;
  for (int i = 0; i < SYM; i++) {
    for (int j = 0; j < SYM; j++) {
      *C_ij = 0.;
      if (i == j) *C_ij += 2. * mu;
      if ((i < DIM) && (j < DIM)) *C_ij += lambda;
      ++C_ij;
    }
  }
  hooke->free = hooke_free;
  hooke->update = hooke_update;
  return hooke;
}

typedef struct HalmDragon1998 {
  struct Material;
  double lambda;
  double mu;
  double alpha; /* Should be > 0 (opposite convention to Xianda's paper. */
  double beta;
  double k0_sqrt2;
  double k1_sqrt2;
} HalmDragon1998;

void halm_dragon_1998_update(HalmDragon1998* mat, double* delta_eps,
                             double* eps1, double* omega1, double* sig2,
                             double* omega2, double* C2) {
  double eps2[SYM];
  for (size_t i = 0; i < SYM; i++) eps2[i] = eps1[i] + delta_eps[i];

  double tr_eps2 = 0.;
  for (size_t i = 0; i < DIM; i++) tr_eps2 += eps2[i];

  double two_alpha_tr_eps2 = 2. * mat->alpha * tr_eps2;
  double four_beta = 4. * mat->beta;
  double H_eps2[SYM];
  for (size_t i = 0; i < SYM; i++) H_eps2[i] = four_beta * eps2[i];
  for (size_t i = 0; i < DIM; i++) H_eps2[i] += two_alpha_tr_eps2;

  double eps2_H_eps2 = 0.;
  for (size_t i = 0; i < SYM; i++) eps2_H_eps2 += eps2[i] * H_eps2[i];

  double f_tr = 0.5 * eps2_H_eps2 - (mat->k0_sqrt2 + mat->k1_sqrt2 * omega1[0]);

  double delta_omega = f_tr > 0. ? f_tr / mat->k1_sqrt2 : 0.;
  omega2[0] = omega1[0] + delta_omega;

  double lambda_tr_eps2 = mat->lambda * tr_eps2;
  double two_mu = 2. * mat->mu;
  for (size_t i = 0; i < SYM; i++)
    sig2[i] = two_mu * eps2[i] - omega2[0] * H_eps2[i];
  for (size_t i = 0; i < DIM; i++) sig2[i] += lambda_tr_eps2;

  if (C2 != NULL) {
    double lambda_sec = mat->lambda - 2. * omega2[0] * mat->alpha;
    double two_mu_sec = 2. * (mat->mu - 2. * omega2[0] * mat->beta);
    double* C2_ij = C2;
    if (f_tr > 0) {
      for (size_t i = 0; i < SYM; i++) {
        double aux = -H_eps2[i] / mat->k1_sqrt2;
        for (size_t j = 0; j < SYM; j++) {
          *C2_ij = aux * H_eps2[j];
          if (i == j) *C2_ij += two_mu_sec;
          if ((i < DIM) && (j < DIM)) *C2_ij += lambda_sec;
          ++C2_ij;
        }
      }
    } else {
      for (size_t i = 0; i < SYM; i++) {
        for (size_t j = 0; j < SYM; j++) {
          *C2_ij = 0;
          if (i == j) *C2_ij += two_mu_sec;
          if ((i < DIM) && (j < DIM)) *C2_ij += lambda_sec;
          ++C2_ij;
        }
      }
    }
  }
}

HalmDragon1998* halm_dragon_1998_new(double lambda, double mu, double alpha,
                                     double beta, double k0, double k1) {
  HalmDragon1998* mat = malloc(sizeof(HalmDragon1998));
  mat->lambda = lambda;
  mat->mu = mu;
  mat->alpha = alpha;
  mat->beta = beta;
  mat->k0_sqrt2 = k0 * M_SQRT2;
  mat->k1_sqrt2 = k1 * M_SQRT2;
  mat->free = free;
  mat->update = halm_dragon_1998_update;
  return mat;
}

HalmDragon1998* halm_dragon_1998_new_default() {
  double kappa = 60700.;
  double mu = 31300.;
  double alpha = 16000.;
  double beta = 31000.;
  double k0 = 0.11;
  double k1 = 2.2;
  return halm_dragon_1998_new(kappa - 2 * mu / DIM, mu, alpha, beta, k0, k1);
}

void hd98_test_hydrostatic_strain() {
  FILE* f = fopen("test.csv", "w");
  HalmDragon1998* mat = halm_dragon_1998_new_default();
  double eps_max = 2.5e-3;
  size_t num_steps = 9;

  double delta_eps[SYM];
  for (size_t j = 0; j < DIM; j++)
    delta_eps[j] = eps_max / (double)(num_steps - 1);

  double eps[] = {0., 0., 0., 0., 0., 0.};
  double omega = 0.;
  double sig[SYM];
  double sig_exp, omega_exp, lambda_exp, mu_exp;
  for (size_t i = 1; i < num_steps; i++) {
    double omega_new;
    mat->update(mat, delta_eps, eps, &omega, sig, &omega_new, NULL);
    for (size_t j = 0; j < SYM; j++) eps[j] += delta_eps[j];
    omega = omega_new;

    /* Expected values */
    omega_exp = (3. * (3. * mat->alpha + 2. * mat->beta) * eps[0] * eps[0] -
                 mat->k0_sqrt2) /
                mat->k1_sqrt2;
    if (omega_exp < 0) omega_exp = 0.;
    sig_exp = (3. * (mat->lambda - 2. * omega_exp * mat->alpha) +
               2 * (mat->mu - 2. * omega_exp * mat->beta)) *
              eps[0];

    fprintf(f, "%g,%g,%g,%g,%g\n", eps[0], sig[0], omega, sig_exp, omega_exp);
    lambda_exp = mat->lambda - 2. * omega_exp * mat->alpha;
    mu_exp = mat->mu - 2. * omega_exp * mat->beta;
  }

  /* Unloading */
  for (size_t j = 0; j < DIM; j++) delta_eps[j] = -delta_eps[j];
  for (size_t i = 1; i < num_steps; i++) {
    double omega_new;
    mat->update(mat, delta_eps, eps, &omega, sig, &omega_new, NULL);
    for (size_t j = 0; j < SYM; j++) eps[j] += delta_eps[j];
    omega = omega_new;

    /* Expected values */
    sig_exp = (3. * lambda_exp + 2 * mu_exp) * eps[0];

    fprintf(f, "%g,%g,%g,%g,%g\n", eps[0], sig[0], omega, sig_exp, omega_exp);
  }

  mat->free(mat);
  fclose(f);
}

void global_update(size_t n, double* delta_eps, double* eps1, double* omega1,
                   uint8_t* phase, Material** mat, double* sig2, double* omega2,
                   double* C2) {
  double* delta_eps_i = delta_eps;
  double* eps1_i = eps1;
  double* omega1_i = omega1;
  Material* mat_i;
  double* sig2_i = sig2;
  double* omega2_i = omega2;
  double* C2_i = C2;

  for (size_t i = 0; i < n; i++) {
    mat_i = mat[phase[i]];
    mat_i->update(mat_i, delta_eps_i, eps1_i, omega1_i, sig2_i, omega2_i, C2_i);

    delta_eps_i += SYM;
    eps1_i += SYM;
    ++omega1_i; /* This assumes that there is only one internal variable.
                 */
    sig2_i += SYM;
    ++omega2_i;
    C2_i += SYM * SYM;
  }
}

void test_global_update() {
  Material** mat = malloc(2 * sizeof(Material*));
  mat[0] = halm_dragon_1998_new_default();
  double kappa1 = 76700.;
  double mu1 = 41600.;
  mat[1] = hooke_new(kappa1 - 2 * mu1 / (double)DIM, mu1);

  size_t n = 10;
  uint8_t* phase = malloc(n * sizeof(size_t));
  for (size_t i = 0; i < n; i++) phase[i] = i % 2;

  double* delta_eps = calloc(n * SYM, sizeof(double));
  double* eps1 = calloc(n * SYM, sizeof(double));
  double* omega1 = calloc(n, sizeof(double));
  double* sig2 = calloc(n * SYM, sizeof(double));
  double* omega2 = calloc(n, sizeof(double));
  double* C2 = calloc(n * SYM * SYM, sizeof(double));

  for (size_t i = 0; i < n; i++) {
    omega1[i] = 0.;
    omega2[i] = 0.;
  }

  global_update(n, delta_eps, eps1, omega1, phase, mat, sig2, omega2, C2);
  double* C2_ijk = C2;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < SYM; j++) {
      for (size_t k = 0; k < SYM; k++) {
        printf("%g\t", *C2_ijk);
        ++C2_ijk;
      }
      printf("\n");
    }
  }

  free(C2);
  free(omega2);
  free(sig2);
  free(omega1);
  free(eps1);
  free(delta_eps);
  free(phase);
  free(mat);
}

int main() {
  test_global_update();
  return 0;
}
