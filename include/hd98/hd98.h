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

typedef struct MaterialType_ MaterialType;
typedef struct Material_ Material;

typedef void material_free_t(Material *);
typedef void material_update_t
    (Material const *, double const *, double const *, double const *,
     double *, double *, double *);

struct MaterialType_ {
  char *name;
  material_free_t *free;
  material_update_t *update;
};

struct Material_ {
  MaterialType const *type;
  void *data;
};

DllExport Material *hd98_hooke_new(double, double);

DllExport Material *hd98_halm_dragon_1998_new(double, double, double, double, double,
                                              double);

DllExport Material *hd98_halm_dragon_1998_new_default();

DllExport void hd98_global_update(size_t n,
                                  double const *delta_eps,
                                  double const *eps1,
                                  double const *omega1,
                                  uint8_t const *phase,
                                  Material **mat,
                                  double *sig2,
                                  double *omega2,
                                  double *C2);

#endif
