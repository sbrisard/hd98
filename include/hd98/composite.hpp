#pragma once

#include <hd98/hd98.hpp>
#include <hd98/hooke.hpp>
#include <hd98/halm_dragon_1998.hpp>

namespace hd98 {
class Composite {
  // This is ugly, but will do for the time being.
  // The truth is that I don't know how to pass vectors of Material to python

 public:
  Hooke const hooke;
  HalmDragon1998 const halm_dragon_1998;

  Composite(Hooke const &hooke, HalmDragon1998 const &halm_dragon_1998)
      : hooke{hooke}, halm_dragon_1998{halm_dragon_1998} {}

  void update(std::span<size_t> phase, double const *delta_eps,
              double const *eps1, double const *omega1, double *sig2,
              double *omega2, double *C2) {
    double const *delta_eps_i = delta_eps;
    double const *eps1_i = eps1;
    double const *omega1_i = omega1;
    double *sig2_i = sig2;
    double *omega2_i = omega2;
    double *C2_i = C2;

    for (auto phase_i : phase) {
      if (phase_i == 0) {
        hooke.update(delta_eps_i, eps1_i, omega1_i, sig2_i, omega2_i, C2_i);
      } else {
        halm_dragon_1998.update(delta_eps_i, eps1_i, omega1_i, sig2_i, omega2_i,
                                C2_i);
      }
      delta_eps_i += sym;
      eps1_i += sym;
      omega1_i += 1;
      sig2_i += sym;
      omega2_i += 1;
      C2_i += sym * sym;
    }
  }
};
}