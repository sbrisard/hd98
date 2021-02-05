#ifndef __TEST_HD98_H_20205717185754__
#define __TEST_HD98_H_20205717185754__

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

void assert_true(bool predicate);
void assert_equal(double act, double exp, double rtol, double atol);
void assert_array_equal(size_t size, double const *actual,
                        double const *expected, double rtol, double atol);
#endif
