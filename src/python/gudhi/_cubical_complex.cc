/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <python_interfaces/Cubical_complex_interface.h>
// #include <python_interfaces/Persistent_cohomology_interface.h>

namespace nb = nanobind;

using CC = Gudhi::cubical_complex::Cubical_complex_interface;
// using CPers = Gudhi::Persistent_cohomology_interface<CC>;

using PCC = Gudhi::cubical_complex::Periodic_cubical_complex_interface;
// using PCPers = Gudhi::Persistent_cohomology_interface<PCC>;

NB_MODULE(_cubical_complex_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<CC>(m, "_Bitmap_cubical_complex_interface")
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, bool>(),
           nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const std::string&>(), nb::call_guard<nb::gil_scoped_release>())
      .def("num_simplices", &CC::num_simplices, nb::call_guard<nb::gil_scoped_release>(), R"doc(
This function returns the number of all cubes in the complex.

:returns:  int -- the number of all cubes in the complex.
           )doc")
      .def("dimension",
           nb::overload_cast<>(&CC::dimension, nb::const_),
           nb::call_guard<nb::gil_scoped_release>(),
           R"doc(
This function returns the dimension of the complex.

:returns:  int -- the complex dimension.
           )doc")
      .def("_shape", &CC::shape)
      .def("_get_numpy_array", &CC::get_numpy_array, nb::rv_policy::reference_internal);

  nb::class_<PCC>(m, "_Periodic_cubical_complex_interface")
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, const std::vector<bool>&, bool>(),
           nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const std::string&>(), nb::call_guard<nb::gil_scoped_release>())
      .def("num_simplices", &PCC::num_simplices, nb::call_guard<nb::gil_scoped_release>(), R"doc(
This function returns the number of all cubes in the complex.

:returns:  int -- the number of all cubes in the complex.
           )doc")
      .def("dimension",
           nb::overload_cast<>(&PCC::dimension, nb::const_),
           nb::call_guard<nb::gil_scoped_release>(),
           R"doc(
This function returns the dimension of the complex.

:returns:  int -- the complex dimension.
           )doc")
      .def("_shape", &PCC::shape)
      .def("_periodicities", &PCC::periodicities)
      .def("_get_numpy_array", &PCC::get_numpy_array, nb::rv_policy::reference_internal);
}
