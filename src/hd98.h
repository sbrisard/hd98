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

typedef void material_free_t(MaterialType *);
typedef void material_update_t(MaterialType *, double *, double *, double *,
                               double *, double *, double *);

struct MaterialType_ {
  material_free_t *free;
  material_update_t *update;
};

typedef struct Material_ {
  MaterialType const *type;
  void *data;
} Material;

DllExport Material *hooke_new(double, double);

DllExport Material *halm_dragon_1998_new(double, double, double, double,
                                               double, double);

DllExport Material *halm_dragon_1998_new_default();

DllExport void global_update(size_t, double *, double *, double *, uint8_t *,
                             MaterialType **, double *, double *, double *);

#endif
