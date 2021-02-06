#pragma once

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

#include <stdint.h>

constexpr size_t HD98_DIM = 3;
constexpr size_t HD98_SYM = 6;

constexpr int HD98_TANGENT_STIFFNESS = 0;
constexpr int HD98_SECANT_STIFFNESS = 1;

namespace hd98 {
typedef struct HD98_MaterialType_ HD98_MaterialType;
typedef struct HD98_Material_ HD98_Material;

typedef void hd98_material_free_t(HD98_Material *mat);
typedef void hd98_material_current_state_t(HD98_Material const *mat,
                                           double const *eps, double const *iv,
                                           double *sig);
typedef void hd98_material_update_t(HD98_Material const *mat,
                                    double const *delta_eps, double const *eps1,
                                    double const *iv1, double *sig2,
                                    double *iv2, double *C2);

struct HD98_MaterialType_ {
  char name[64]; /* Modify core.py if this length is altered. */
  size_t niv;    /* Number of internal variables. */
  hd98_material_free_t *free;
  hd98_material_current_state_t *current_state;
  hd98_material_update_t *update;
};

struct HD98_Material_ {
  HD98_MaterialType const *type;
  void *data;
};

DllExport void hd98_global_update(size_t n, size_t const *phase,
                                  HD98_Material const **mat,
                                  double const *delta_eps, double const *eps1,
                                  double const *iv1, double *sig2, double *iv2,
                                  double *C2);
}  // namespace hd98