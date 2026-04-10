/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <gudhi/vineyard_helper.h>
#include <python_interfaces/Persistent_cohomology_interface.h>
#include <python_interfaces/complex_pcoh_interfaces.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/Cubical_complex_interface.h>

namespace nb = nanobind;

namespace Gudhi {

template <class FilteredComplex>
inline nb::tuple build_python_boundary_matrix_from_complex(FilteredComplex &complex, int maxDimension) {
  std::vector<std::vector<typename FilteredComplex::Simplex_key>> boundaries;
  std::vector<int> dimensions;
  std::vector<typename FilteredComplex::Filtration_value> filtrationValues;

  Gudhi::vineyard::build_boundary_matrix_from_complex(complex, boundaries, dimensions, filtrationValues, maxDimension);

  nb::list b;
  for (auto &boundary : boundaries) {
    b.append(_wrap_as_numpy_array(std::move(boundary), boundary.size()));
  }
  return nb::make_tuple(std::move(b), _wrap_as_numpy_array(std::move(dimensions), dimensions.size()),
                        _wrap_as_numpy_array(std::move(filtrationValues), filtrationValues.size()));
}

}  // namespace Gudhi

using gcpi = Gudhi::Complex_pcoh_interface;
using gsti = Gudhi::Simplex_tree_interface;
using gcci = Gudhi::cubical_complex::Cubical_complex_interface;
using gpcci = Gudhi::cubical_complex::Periodic_cubical_complex_interface;

NB_MODULE(_build_helpers_ext, m) {
  m.attr("__license__") = "MIT";

  m.def("_build_boundary_matrix_from_complex", &Gudhi::build_python_boundary_matrix_from_complex<gsti>);
  m.def("_build_boundary_matrix_from_complex", &Gudhi::build_python_boundary_matrix_from_complex<gcci>);
  m.def("_build_boundary_matrix_from_complex", &Gudhi::build_python_boundary_matrix_from_complex<gpcci>);
}
