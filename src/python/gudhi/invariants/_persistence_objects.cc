/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>  // std::size_t
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

nb::ndarray<nb::numpy, int> get_betti_numbers(const nb::list& finites, const nb::list& infinites, double min,
                                                 double max) {
  std::vector<int> bettiNumbers(finites.size());

  int currDim = 0;
  for (nb::handle h : finites) {
    auto finites_view = nb::cast<nb::ndarray<const double, nb::ndim<2>>>(h).view();
    for (std::size_t i = 0; i < finites_view.shape(0); ++i) {
      if (finites_view(i, 0) <= min && finites_view(i, 1) > max) ++bettiNumbers[currDim];
    }
    ++currDim;
  }

  currDim = 0;
  for (nb::handle h : infinites) {
    auto infinites_view = nb::cast<nb::ndarray<const double, nb::ndim<1>>>(h).view();
    for (std::size_t i = 0; i < infinites_view.shape(0); ++i) {
      if (infinites_view(i) <= min) ++bettiNumbers[currDim];
    }
    ++currDim;
  }

  return _wrap_as_numpy_array(std::move(bettiNumbers), bettiNumbers.size());
}

NB_MODULE(_persistence_objects_ext, m) {
  m.attr("__license__") = "MIT";
  m.def("_get_betti_numbers", &get_betti_numbers);
}
