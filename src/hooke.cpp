#include <ostream>

#include "hd98/hooke.hpp"

namespace hd98{
std::ostream &operator<<(std::ostream &os, const Hooke mat) {
  return os << mat.repr();
}
}