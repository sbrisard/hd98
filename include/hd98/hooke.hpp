#ifndef __HOOKE_H_20203317183340__
#define __HOOKE_H_20203317183340__

#include "hd98/hd98.hpp"

typedef struct HD98_HookeData_ {
  double lambda;
  double mu;
  double *C;
} HD98_HookeData;

DllExport HD98_Material *hd98_hooke_new(double lambda, double mu);

#endif
