#pragma once
#include <ostream>
#include <string>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

namespace hd98 {
constexpr size_t dim = 3;
constexpr size_t sym = 6;

constexpr int tangent_stiffness = 0;
constexpr int secant_stiffness = 1;

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

class Material {
 public:
  [[nodiscard]] virtual std::string repr() const = 0;
  virtual void current_state(double const *eps, double const *omega,
                             double *sig) const = 0;
  virtual void update(double const *delta_eps, double const *eps1,
                      double const *omega1, double *sig2, double *omega2,
                      double *C2) const = 0;
};

std::ostream &operator<<(std::ostream &os, const Material &mat);

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