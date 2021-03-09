#include <array>
#include <iostream>
#include "hd98/hd98.hpp"
#include "hd98/hooke.hpp"

int main() {
  double const lambda = 1.2;
  double const mu = 1.0;
  hd98::Hooke hooke{lambda, mu};
  std::cout << hooke << std::endl;
  std::array<double, hd98::sym> eps{-1.2, 3.4, -5.6, 7.8, -9.1, 10.11};
  std::array<double, hd98::sym> act;
  hooke.current_state(eps.data(), nullptr, act.data());

  std::array<double, hd98::sym> exp;
  double lambda_tr_eps = lambda * (eps[0] + eps[1] + eps[2]);
  for (int i = 0; i < hd98::sym; i++) {
    exp[i] = 2 * hooke.mu * eps[i];
  }
  for (int i = 0; i < hd98::dim; i++) {
    exp[i] += lambda_tr_eps;
  }
  for (int i = 0; i < hd98::sym; i++) {
    std::cout << "sig[" << i << "] = " << exp[i] << " (expected) = " << act[i]
              << " (actual)" << std::endl;
  }
}
