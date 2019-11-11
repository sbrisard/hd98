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

/* General material model. */

typedef struct Material_ Material;

typedef void material_free_t(Material *);
typedef void material_update_t(Material *, double *, double *, double *,
                               double *, double *, double *);

struct Material_ {
  material_free_t *free;
  material_update_t *update;
};

/* Hooke's model. */

typedef struct Hooke {
  struct Material_;
  double lambda;
  double mu;
  double *C;
} Hooke;

DllExport Hooke *hooke_new(double, double);

/* Model of Halm and Dragon (1998). */

typedef struct HalmDragon1998 {
  struct Material_;
  double lambda;
  double mu;
  double alpha; /* Should be > 0 (opposite convention to Xianda's paper. */
  double beta;
  double k0_sqrt2;
  double k1_sqrt2;
} HalmDragon1998;

DllExport HalmDragon1998 *halm_dragon_1998_new(double, double, double, double,
                                               double, double);

DllExport HalmDragon1998 *halm_dragon_1998_new_default();

DllExport void global_update(size_t, double *, double *, double *, uint8_t *,
                             Material **, double *, double *, double *);

#endif
