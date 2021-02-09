#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "hd98/halm_dragon_1998.hpp"
#include "hd98/hd98.hpp"
#include "hd98/hooke.hpp"

namespace py = pybind11;

using array1d = py::array_t<double>;

PYBIND11_MODULE(pyhd98, m) {
  m.doc() = "Python bindings to the hd98 library";
  m.attr("__author__") = pybind11::cast(__HD98_AUTHOR__);
  m.attr("__version__") = pybind11::cast(__HD98_VERSION__);

  py::class_<hd98::Hooke>(m, "Hooke")
      .def(py::init<double, double>())
      .def("__repr__", &(hd98::Hooke::repr))
      .def_readonly("lambda_", &hd98::Hooke::lambda)
      .def_readonly("mu", &hd98::Hooke::mu)
      .def("current_state",
           [](hd98::Hooke& self, array1d eps, array1d unused, array1d sig) {
             return self.current_state(eps.data(), nullptr, sig.mutable_data());
           })
      .def("update", [](hd98::Hooke& self, array1d delta_eps, array1d eps1,
                        array1d unused, array1d sig2, array1d unused2,
                        array1d C2) {
        return self.update(delta_eps.data(), eps1.data(), nullptr,
                           sig2.mutable_data(), nullptr, C2.mutable_data());
      });

  py::class_<hd98::HalmDragon1998>(m, "HalmDragon1998")
      .def(py::init<double, double, double, double, double, double, int>())
      .def("__repr__", &(hd98::HalmDragon1998::repr))
      .def_readonly("lambda_", &hd98::HalmDragon1998::lambda)
      .def_readonly("mu", &hd98::HalmDragon1998::mu)
      .def_readonly("alpha", &hd98::HalmDragon1998::alpha)
      .def_readonly("beta", &hd98::HalmDragon1998::beta)
      .def_readonly("k0_sqrt2", &hd98::HalmDragon1998::k0_sqrt2)
      .def_readonly("k1_sqrt2", &hd98::HalmDragon1998::k1_sqrt2)
      .def_readonly("stiffness_type", &hd98::HalmDragon1998::stiffness_type)
      .def("current_state",
           [](hd98::HalmDragon1998& self, array1d eps, array1d omega,
              array1d sig) {
             return self.current_state(eps.data(), omega.data(), sig.mutable_data());
           })
      .def("update",
           [](hd98::HalmDragon1998& self, array1d delta_eps, array1d eps1,
              array1d omega1, array1d sig2, array1d omega2, array1d C2) {
             return self.update(delta_eps.data(), eps1.data(), omega1.data(),
                                sig2.mutable_data(), omega2.mutable_data(),
                                C2.mutable_data());
           });
}
