#ifndef __HALM_DRAGON_1998_H_20204117184143__
#define __HALM_DRAGON_1998_H_20204117184143__

#include "hd98/hd98.h"

typedef struct HD98_HalmDragon1998Data_ {
  double lambda;
  double mu;
  double alpha; /* Should be > 0 (opposite convention to Xianda's paper. */
  double beta;
  double k0_sqrt2;
  double k1_sqrt2;
  int stiffness_type;
} HD98_HalmDragon1998Data;

HD98_MaterialType const HD98_HalmDragon1998;

DllExport HD98_Material *hd98_halm_dragon_1998_new(double lambda, double mu,
                                                   double alpha, double beta,
                                                   double k0, double k1,
                                                   int stiffness_type);

DllExport HD98_Material *hd98_halm_dragon_1998_new_default();

#endif
