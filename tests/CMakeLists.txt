add_executable(test_hd98
        test_hd98.c
        test_hooke.c
        test_halm_dragon_1998.c
        )

target_include_directories(test_hd98 PUBLIC ../include)

find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(test_hd98 ${MATH_LIBRARY})
endif ()

target_link_libraries(test_hd98 PRIVATE hd98)

#find_package(Catch2 REQUIRED)
#target_link_libraries(test_hd98 PRIVATE hd98 Catch2::Catch2)

add_test(NAME test_hd98 COMMAND test_hd98)
