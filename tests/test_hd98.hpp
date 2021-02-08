#ifndef __TEST_HD98_H_20205717185754__
#define __TEST_HD98_H_20205717185754__

#include <array>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "hd98/hd98.hpp"

using Tensor2 = std::array<double, hd98::sym>;
using Tensor4 = std::array<double, hd98::sym * hd98::sym>;

void assert_equal(double act, double exp, double rtol, double atol);
void assert_array_equal(size_t size, double const *actual,
                        double const *expected, double rtol, double atol);
#endif
