#include "hd98/hooke.hpp"

namespace hd98 {
std::array<double, HD98_SYM * HD98_SYM> stiffness_matrix(double lambda,
                                                         double mu) {
  std::array<double, HD98_SYM * HD98_SYM> C{};
  for (size_t i = 0, ij = 0; i < HD98_SYM; i++) {
    for (size_t j = 0; j < HD98_SYM; j++, ij++) {
      C[ij] = 0.0;
      if (i == j) C[ij] += 2. * mu;
      if ((i < HD98_DIM) && (j < HD98_DIM)) C[ij] += lambda;
    }
  }
  return C;
}
}  // namespace hd98