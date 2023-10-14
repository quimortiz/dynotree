
#include "dynotree/KDTree.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

template <typename T>
void declare_tree(py::module &m, const std::string &name) {
  py::class_<typename T::DistanceId>(m, (name + "dp").c_str())
      .def(py::init<>())
      .def_readonly("distance", &T::DistanceId::distance)
      .def_readonly("id", &T::DistanceId::id);

  py::class_<T>(m, name.c_str())
      .def(py::init<int>())
      .def("addPoint", &T::addPoint)
      .def("search", &T::search)
      .def("searchKnn", &T::searchKnn)
      .def("getStateSpace", &T::getStateSpace);
  // .def("interpolate", &T::interpolate);
  //
}

template <typename T>
void declare_treex(py::module &m, const std::string &name) {
  // const std::string name = "TreeX";
  py::class_<typename T::DistanceId>(m, (name + "dp").c_str())
      .def(py::init<>())
      .def_readonly("distance", &T::DistanceId::distance)
      .def_readonly("id", &T::DistanceId::id);

  py::class_<T>(m, name.c_str())
      .def(py::init<int, const std::vector<std::string>>())
      .def("addPoint", &T::addPoint)            // add point
      .def("search", &T::search)                // search
      .def("searchKnn", &T::searchKnn)          // search
      .def("getStateSpace", &T::getStateSpace); //
  //
  //
  //
  // .def("interpolate", &T::interpolate);
  //
}

template <typename T>
void declare_state_space_x(py::module &m, const std::string &name) {

  py::class_<T>(m, name.c_str())
      .def(py::init<const std::vector<std::string>>())
      // .def(py::init<>())
      .def("interpolate", &T::interpolate)       // add point
      .def("set_bounds", &T::set_bounds)         // search
      .def("sample_uniform", &T::sample_uniform) // search
      .def("distance", &T::distance)
      .def("distance_to_rectangle", &T::distance_to_rectangle);

  // search
  //
}

template <typename T>
void declare_state_space(py::module &m, const std::string &name) {

  py::class_<T>(m, name.c_str())
      // .def(py::init<int, const std::vector<std::string>>())
      .def(py::init<>())
      .def("interpolate", &T::interpolate)       // add point
      .def("set_bounds", &T::set_bounds)         // search
      .def("sample_uniform", &T::sample_uniform) // search
      .def("distance", &T::distance)
      .def("distance_to_rectangle", &T::distance_to_rectangle);

  // search
  //
}

PYBIND11_MODULE(pydynotree, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  using R2 = dynotree::Rn<double, 2>;
  using R3 = dynotree::Rn<double, 3>;
  using R4 = dynotree::Rn<double, 4>;
  using R5 = dynotree::Rn<double, 5>;
  using R6 = dynotree::Rn<double, 6>;
  using R7 = dynotree::Rn<double, 7>;
  using RX = dynotree::Rn<double, -1>;

  using SO2 = dynotree::SO2<double>;
  using SO3 = dynotree::SO3<double>;

  using R2SO2 = dynotree::R2SO2<double>;
  using R3SO3 = dynotree::R3SO3<double>;

  using Combined = dynotree::Combined<double>;

  declare_state_space<R2>(m, "R2");
  declare_state_space<R3>(m, "R3");
  declare_state_space<R4>(m, "R4");
  declare_state_space<R5>(m, "R5");
  declare_state_space<R6>(m, "R6");
  declare_state_space<R7>(m, "R7");
  declare_state_space<RX>(m, "RX");

  declare_state_space<R2SO2>(m, "R2SO2");
  declare_state_space<R3SO3>(m, "R3SO3");

  declare_state_space<SO3>(m, "SO3");
  declare_state_space<SO2>(m, "SO2");

  declare_state_space_x<Combined>(m, "SpaceX");

  const int bucket_size = 32;
  using TreeRX = dynotree::KDTree<int, -1, bucket_size, double, RX>;
  using TreeR2 = dynotree::KDTree<int, 2, bucket_size, double, R2>;
  using TreeR4 = dynotree::KDTree<int, 4, bucket_size, double, R4>;
  using TreeR7 = dynotree::KDTree<int, 7, bucket_size, double, R7>;
  using TreeR2SO2 = dynotree::KDTree<int, 3, bucket_size, double, R2SO2>;
  using TreeSO3 = dynotree::KDTree<int, 4, bucket_size, double, SO3>;
  using TreeSO2 = dynotree::KDTree<int, 1, bucket_size, double, SO2>;
  using TreeR3SO3 = dynotree::KDTree<int, 7, bucket_size, double, R3SO3>;
  using TreeX = dynotree::KDTree<int, -1, bucket_size, double, Combined>;

  declare_tree<TreeRX>(m, "TreeRX");
  declare_tree<TreeR2>(m, "TreeR2");
  declare_tree<TreeR4>(m, "TreeR4");
  declare_tree<TreeR7>(m, "TreeR7");
  declare_tree<TreeR2SO2>(m, "TreeR2SO2");
  declare_tree<TreeSO3>(m, "TreeSO3");
  declare_tree<TreeSO2>(m, "TreeSO2");
  declare_tree<TreeR3SO3>(m, "TreeR3SO3");
  declare_treex<TreeX>(m, "TreeX");

  m.def("srand", [](int seed) { srand(seed); });
}
