#ifndef __HD98_H_201910301551__
#define __HD98_H_201910301551__

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

#include <stdint.h>

#define HD98_DIM 3
#define HD98_SYM 6

typedef struct HD98_MaterialType_ HD98_MaterialType;
typedef struct HD98_Material_ HD98_Material;

typedef void hd98_material_free_t(HD98_Material *);
typedef void hd98_material_update_t(HD98_Material const *, double const *,
                                    double const *, double const *, double *,
                                    double *, double *);

struct HD98_MaterialType_ {
  char *name;
  hd98_material_free_t *free;
  hd98_material_update_t *update;
};

struct HD98_Material_ {
  HD98_MaterialType const *type;
  void *data;
};

DllExport HD98_Material *hd98_hooke_new(double, double);

DllExport HD98_Material *hd98_halm_dragon_1998_new(double, double, double,
                                                   double, double, double);

DllExport HD98_Material *hd98_halm_dragon_1998_new_default();

DllExport void hd98_global_update(size_t n, double const *delta_eps,
                                  double const *eps1, double const *omega1,
                                  size_t const *phase, HD98_Material **mat,
                                  double *sig2, double *omega2, double *C2);

DllExport int hd98_solve_polarization_plus(HD98_Material *mat, double lambda0,
                                           double mu0, double const *delta_tau,
                                           double const *eps1,
                                           double const *omega1,
                                           double *delta_eps);

DllExport int hd98_solve_polarizations_plus(size_t n, size_t const *phase,
                                            HD98_Material **mat, double lambda0,
                                            double mu0, double const *delta_tau,
                                            double const *eps1,
                                            double const *omega1,
                                            double *delta_eps);
#endif
