find_package(PkgConfig REQUIRED)
pkg_check_modules(deps REQUIRED IMPORTED_TARGET glib-2.0)

add_executable(test_hd98 test_hd98.c)
target_link_libraries(test_hd98 PRIVATE hd98 PkgConfig::deps)
add_test(NAME test_hd98 COMMAND test_hd98)
