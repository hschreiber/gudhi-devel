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

#include <python_interfaces/Persistent_cohomology_interface.h>
#include <python_interfaces/complex_pcoh_interfaces.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/Cubical_complex_interface.h>

namespace nb = nanobind;

using gcpi = Gudhi::Complex_pcoh_interface;
using gdtpi_gcpi = Gudhi::Dimension_threshold_pcoh_interface<gcpi>;
using gpci_gcpi = Gudhi::Persistent_cohomology_interface<gcpi>;
using gpci_gdtpi_gcpi = Gudhi::Persistent_cohomology_interface<gdtpi_gcpi>;

using gsti = Gudhi::Simplex_tree_interface;
using gdtpi_gsti = Gudhi::Dimension_threshold_pcoh_interface<gsti>;
using gpci_gsti = Gudhi::Persistent_cohomology_interface<gsti>;
using gpci_gdtpi_gsti = Gudhi::Persistent_cohomology_interface<gdtpi_gsti>;

using gcci = Gudhi::cubical_complex::Cubical_complex_interface;
using gdtpi_gcci = Gudhi::Dimension_threshold_pcoh_interface<gcci>;
using gpci_gcci = Gudhi::Persistent_cohomology_interface<gcci>;
using gpci_gdtpi_gcci = Gudhi::Persistent_cohomology_interface<gdtpi_gcci>;

using gpcci = Gudhi::cubical_complex::Periodic_cubical_complex_interface;
using gdtpi_gpcci = Gudhi::Dimension_threshold_pcoh_interface<gpcci>;
using gpci_gpcci = Gudhi::Persistent_cohomology_interface<gpcci>;
using gpci_gdtpi_gpcci = Gudhi::Persistent_cohomology_interface<gdtpi_gpcci>;

NB_MODULE(_pers_cohomology_ext, m) {
  m.attr("__license__") = "MIT";
  nb::class_<gdtpi_gcpi>(m, "_Complex_dim_overlay").def(nb::init<gcpi &, int>());
  nb::class_<gdtpi_gsti>(m, "_St_dim_overlay").def(nb::init<gsti &, int>(), nb::call_guard<nb::gil_scoped_release>());
  nb::class_<gdtpi_gcci>(m, "_Cc_dim_overlay").def(nb::init<gcci &, int>(), nb::call_guard<nb::gil_scoped_release>());
  nb::class_<gdtpi_gpcci>(m, "_Pcc_dim_overlay")
      .def(nb::init<gpcci &, int>(), nb::call_guard<nb::gil_scoped_release>());

  nb::class_<gcpi>(m, "_complex_coho_overlay").def(nb::init<nanobind::object &>());
  nb::class_<gpci_gcpi>(m, "_Complex_persistence_interface")
      .def(nb::init<gcpi &, bool>())
      .def("_compute_persistence", &gpci_gcpi::compute_persistence)
      .def("_get_intervals", &gpci_gcpi::build_intervals)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gcpi::build_simplicial_intervals_as_vertices);
  nb::class_<gpci_gdtpi_gcpi>(m, "_Complex_max_dim_persistence_interface")
      .def(nb::init<gdtpi_gcpi &, bool>())
      .def("_compute_persistence", &gpci_gdtpi_gcpi::compute_persistence)
      .def("_get_intervals", &gpci_gdtpi_gcpi::build_intervals)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gdtpi_gcpi::build_simplicial_intervals_as_vertices);

  nb::class_<gpci_gsti>(m, "_Simplex_tree_persistence_interface")
      .def(nb::init<gsti &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gsti::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_persistence", &gpci_gsti::get_persistence)
      .def("_get_intervals", &gpci_gsti::build_intervals)
      .def("_betti_numbers", &gpci_gsti::betti_numbers)
      .def("_persistent_betti_numbers", &gpci_gsti::persistent_betti_numbers)
      .def("_intervals_in_dimension", &gpci_gsti::intervals_in_dimension)
      .def("_write_output_diagram", &gpci_gsti::write_output_diagram, nb::call_guard<nb::gil_scoped_release>())
      .def("_persistence_pairs", &gpci_gsti::persistence_pairs)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gsti::build_simplicial_intervals_as_vertices)
      .def("_lower_star_generators", &gpci_gsti::lower_star_generators)
      .def("_flag_generators", &gpci_gsti::flag_generators)
      .def("_compute_extended_persistence_subdiagrams", &gpci_gsti::compute_extended_persistence_subdiagrams);
  nb::class_<gpci_gdtpi_gsti>(m, "_Simplex_tree_max_dim_persistence_interface")
      .def(nb::init<gdtpi_gsti &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gdtpi_gsti::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_intervals", &gpci_gdtpi_gsti::build_intervals)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gdtpi_gsti::build_simplicial_intervals_as_vertices);

  nb::class_<gpci_gcci>(m, "_Cubical_complex_persistence_interface")
      .def(nb::init<gcci &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gcci::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_persistence", &gpci_gcci::get_persistence)
      .def("_get_intervals", &gpci_gcci::build_intervals)
      .def("_cofaces_of_cubical_persistence_pairs", &gpci_gcci::cofaces_of_cubical_persistence_pairs,
           nb::call_guard<nb::gil_scoped_release>())
      .def("_vertices_of_cubical_persistence_pairs", &gpci_gcci::vertices_of_cubical_persistence_pairs,
           nb::call_guard<nb::gil_scoped_release>())
      .def("_betti_numbers", &gpci_gcci::betti_numbers)
      .def("_persistent_betti_numbers", &gpci_gcci::persistent_betti_numbers)
      .def("_intervals_in_dimension", &gpci_gcci::intervals_in_dimension)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gcci::build_simplicial_intervals_as_vertices);
  nb::class_<gpci_gdtpi_gcci>(m, "_Cubical_complex_max_dim_persistence_interface")
      .def(nb::init<gdtpi_gcci &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gdtpi_gcci::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_intervals", &gpci_gdtpi_gcci::build_intervals)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gdtpi_gcci::build_simplicial_intervals_as_vertices);

  nb::class_<gpci_gpcci>(m, "_Periodic_cubical_complex_persistence_interface")
      .def(nb::init<gpcci &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gpcci::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_persistence", &gpci_gpcci::get_persistence)
      .def("_get_intervals", &gpci_gpcci::build_intervals)
      .def("_cofaces_of_cubical_persistence_pairs", &gpci_gpcci::cofaces_of_cubical_persistence_pairs,
           nb::call_guard<nb::gil_scoped_release>())
      .def("_vertices_of_cubical_persistence_pairs", &gpci_gpcci::vertices_of_cubical_persistence_pairs,
           nb::call_guard<nb::gil_scoped_release>())
      .def("_betti_numbers", &gpci_gpcci::betti_numbers)
      .def("_persistent_betti_numbers", &gpci_gpcci::persistent_betti_numbers)
      .def("_intervals_in_dimension", &gpci_gpcci::intervals_in_dimension)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gpcci::build_simplicial_intervals_as_vertices);
  nb::class_<gpci_gdtpi_gpcci>(m, "_Periodic_cubical_complex_max_dim_persistence_interface")
      .def(nb::init<gdtpi_gpcci &, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &gpci_gdtpi_gpcci::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_intervals", &gpci_gdtpi_gpcci::build_intervals)
      .def("_get_simplicial_intervals_as_vertices", &gpci_gdtpi_gpcci::build_simplicial_intervals_as_vertices);
}
