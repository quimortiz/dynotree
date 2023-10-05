
#include "dynkdtree/KDTree.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>



int add(int i, int j) { return i + j; }

// PYBIND11_MODULE(example, m) {
//     declare_foo<int>(m, "Int");
//     declare_foo<double>(m, "Double");
// }

namespace py = pybind11;

template <typename T>
void declare_tree(py::module &m, const std::string &name) {
  py::class_<typename T::DistancePayload>(m, (name + "dp").c_str())
      .def(py::init<>())
      .def_readonly("distance", &T::DistancePayload::distance)
      .def_readonly("payload", &T::DistancePayload::payload);

  py::class_<T>(m, name.c_str())
      .def(py::init<int>())            // constructor
      .def("addPoint", &T::addPoint)   // add point
      .def("search", &T::search)       // search
      .def("searchKnn", &T::searchKnn) // search
      .def("interpolate", &T::interpolate);
  //
}

PYBIND11_MODULE(dynotree, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  m.def("add", &add, "A function that adds two numbers");

  using TreeRX = jk::tree::KDTree<int, -1>;
  using TreeR2 = jk::tree::KDTree<int, 2>;
  using TreeR4 = jk::tree::KDTree<int, 4>;
  using TreeR7 = jk::tree::KDTree<int, 7>;
  using TreeR2SO2 =
      jk::tree::KDTree<int, 3, 23, double, jk::tree::R2SO2<double>>;

  declare_tree<TreeRX>(m, "TreeRX");
  declare_tree<TreeR2>(m, "TreeR2");
  declare_tree<TreeR4>(m, "TreeR4");
  declare_tree<TreeR7>(m, "TreeR7");
  declare_tree<TreeR2SO2>(m, "TreeR2SO2");
}
