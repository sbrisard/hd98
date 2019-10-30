project('hd98', 'c')
math = meson.get_compiler('c').find_library('m', required : false)

libhd98 = library('hd98', ['hd98.c', 'hd98.h'], dependencies: math)
