cmake_minimum_required(VERSION 3.13)
project(hd98 C)

set(CMAKE_C_STANDARD 11)
include(CTest)

find_library(MATH_LIBRARY m)

add_library(hd98 SHARED src/hd98.c include/hd98/hd98.h)
target_include_directories(hd98 PUBLIC include)
if (MATH_LIBRARY)
  target_link_libraries(hd98 ${MATH_LIBRARY})
endif()

add_subdirectory(tests)
