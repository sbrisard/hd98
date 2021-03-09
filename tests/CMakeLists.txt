add_executable(test_hd98
        test_hd98_main.cpp
        test_hd98.cpp
        )

target_include_directories(test_hd98 PUBLIC ../include)

find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(test_hd98 ${MATH_LIBRARY})
endif ()
add_compile_definitions(_USE_MATH_DEFINES)

target_link_libraries(test_hd98 PRIVATE hd98)

find_package(Catch2 REQUIRED)
target_link_libraries(test_hd98 PRIVATE hd98 Catch2::Catch2)

add_test(NAME test_hd98 COMMAND test_hd98)
