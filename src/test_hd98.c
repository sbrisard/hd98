#include <hd98.h>

void test_hd98_proportional_strain() {
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

int main(int argc, char** argv) {
  test_hd98_proportional_strain();
  test_global_update();
}
