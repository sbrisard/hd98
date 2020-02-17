#include <glib.h>

#include "hd98/hooke.h"
#include "test_hd98.h"

void test_material_type() {
  g_assert_cmpstr(HD98_Hooke.name, ==, "Hooke");
  g_assert_cmpuint(HD98_Hooke.num_int_var, ==, 0);
}

void setup_hooke_tests() {
    g_test_add_func("/Hooke/HD98_MaterialType", test_material_type);
}
