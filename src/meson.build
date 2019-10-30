project('hd98', 'c')
math = meson.get_compiler('c').find_library('m', required : false)

hd98_lib = library('hd98', ['hd98.c', 'hd98.h'], dependencies: math)

hd98_test = executable('test_hd98', 'test_hd98.c',
		       dependencies: math,
		       link_with: hd98_lib)
