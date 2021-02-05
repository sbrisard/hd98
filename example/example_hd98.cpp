#include <iostream>
#include "hd98/hd98.hpp"

int main() {
  std::cout << "version = " << hd98::version() << std::endl;
  std::cout << "author = " << hd98::author() << std::endl;
  std::cout << "return_one() = " << hd98::return_one() << std::endl;
  return 0;
}
