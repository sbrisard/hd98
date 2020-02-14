find_library(MATH_LIBRARY m)

add_library(hd98 SHARED hd98.c ../include/hd98/hd98.h)

target_include_directories(hd98 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}../include>
  $<INSTALL_INTERFACE:include>)

target_include_directories(hd98 PUBLIC ../include)
if (MATH_LIBRARY)
  target_link_libraries(hd98 ${MATH_LIBRARY})
endif()

include(FindGSL)
target_link_libraries(hd98 GSL::gsl)

install(TARGETS hd98 LIBRARY DESTINATION lib)
install(DIRECTORY ../include/hd98 DESTINATION include)
