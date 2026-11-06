/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Nested_array_filtration.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/complex_helpers.h>

using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_filtration::Nested_array_filtration;
using Gudhi::multi_persistence::Box;
using Gudhi::multi_persistence::get_bounding_box_of_complex;
using Gudhi::multi_persistence::get_vertices_from_boundaries;
using Gudhi::multi_persistence::get_vertices_from_ordered_boundaries;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::normalize_filtration_values_in_complex;
using Gudhi::multi_persistence::null_vertex;

BOOST_AUTO_TEST_CASE(complex_helpers_get_vertices_ordered) {
  using Idx = std::uint32_t;
  using Vec = std::vector<std::vector<Idx>>;

  // vertices 0,1,2 ; edges 3={0,1}, 4={1,2}, 5={0,2} ; triangle 6
  Vec bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}};
  std::vector<int> dc = {0, 0, 0, 1, 1, 1, 2};

  auto v = get_vertices_from_ordered_boundaries<Idx>(bc, dc);
  Vec expected = {{0}, {1}, {2}, {0, 1}, {1, 2}, {0, 2}, {0, 1, 2}};
  BOOST_CHECK(v == expected);

  // maxDim = 1: triangle is dropped
  auto v1 = get_vertices_from_ordered_boundaries<Idx>(bc, dc, 1);
  BOOST_CHECK_EQUAL(v1.size(), 6);
  BOOST_CHECK(v1 == Vec(expected.begin(), expected.begin() + 6));

  // maxDim = 0: only vertices
  auto v0 = get_vertices_from_ordered_boundaries<Idx>(bc, dc, 0);
  BOOST_CHECK(v0 == Vec(expected.begin(), expected.begin() + 3));

  // maxDim greater than the max dimension: everything
  auto v5 = get_vertices_from_ordered_boundaries<Idx>(bc, dc, 5);
  BOOST_CHECK(v5 == expected);

  // empty input
  Vec bcEmpty;
  std::vector<int> dcEmpty;
  BOOST_CHECK(get_vertices_from_ordered_boundaries<Idx>(bcEmpty, dcEmpty, -1).empty());

  // no dimension 0 cell: every cell is mapped to {null_vertex}
  Vec bc2 = {{}, {0}};
  std::vector<int> dc2 = {1, 2};
  auto v2 = get_vertices_from_ordered_boundaries<Idx>(bc2, dc2);
  BOOST_CHECK_EQUAL(v2.size(), 2);
  BOOST_CHECK_EQUAL(v2[0].size(), 1);
  BOOST_CHECK_EQUAL(v2[0][0], null_vertex<Idx>());
  BOOST_CHECK_EQUAL(v2[1][0], null_vertex<Idx>());
}

BOOST_AUTO_TEST_CASE(complex_helpers_get_vertices_general) {
  using Idx = std::uint32_t;
  using Vec = std::vector<std::vector<Idx>>;

  // cells are not ordered by dimension, with interleaved vertices
  // 0: edge {1,2}, 1: vertex, 2: vertex, 3: triangle, 4: edge {1,3}, 5: vertex, 6: edge {2,3}
  Vec bc = {{1, 2}, {}, {}, {0, 4, 6}, {1, 5}, {}, {2, 5}};
  std::vector<int> dc = {1, 0, 0, 2, 1, 0, 1};
  auto v = get_vertices_from_boundaries<Idx>(bc, dc);
  // vertices are numbered in order of appearance: cell 1 -> 0, cell 2 -> 1, cell 5 -> 2
  Vec expected = {{0, 1}, {0}, {1}, {0, 1, 2}, {0, 2}, {2}, {1, 2}};
  BOOST_CHECK(v == expected);

  // ordered case gives the same answer as the ordered function
  Vec bc2 = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}};
  std::vector<int> dc2 = {0, 0, 0, 1, 1, 1, 2};
  BOOST_CHECK(get_vertices_from_boundaries<Idx>(bc2, dc2) == get_vertices_from_ordered_boundaries<Idx>(bc2, dc2));

  // maxDim: cells above the dimension are not computed (left empty)
  auto v1 = get_vertices_from_boundaries<Idx>(bc2, dc2, 1);
  BOOST_CHECK_EQUAL(v1.size(), 7);
  BOOST_CHECK(v1[6].empty());
  BOOST_CHECK((v1[3] == std::vector<Idx>{0, 1}));

  // empty-boundary cells of dimension > 0 get the null vertex
  Vec bc3 = {{}, {}, {0, 1}};
  std::vector<int> dc3 = {0, 1, 2};
  auto v3 = get_vertices_from_boundaries<Idx>(bc3, dc3);
  BOOST_CHECK_EQUAL(v3[0][0], 0);
  BOOST_CHECK_EQUAL(v3[1][0], null_vertex<Idx>());

  // cycles are detected
  Vec bcCycle = {{}, {2}, {1}};
  std::vector<int> dcCycle = {0, 1, 1};
  BOOST_CHECK_THROW(get_vertices_from_boundaries<Idx>(bcCycle, dcCycle), std::invalid_argument);
}

template <class Fil>
using Complex_t = Multi_parameter_filtered_complex<Fil, std::uint32_t, int>;

using list_of_tested_variants = boost::mpl::list<Multi_parameter_filtration_value<Flat_array_filtration<double>>,
                                                 Multi_parameter_filtration_value<Nested_array_filtration<double>>>;

template <class Fil>
Complex_t<Fil> build_complex() {
  using C = Complex_t<Fil>;
  using ini = std::initializer_list<typename Fil::value_type>;
  typename C::Boundary_container bc = {{}, {}, {}, {0, 1}, {1, 2}};
  typename C::Dimension_container dc = {0, 0, 0, 1, 1};
  typename C::Filtration_value_container fc = {ini{0, 4}, ini{2, 1}, ini{1, 3}, ini{2, 4}, ini{4, 3}};
  return C(bc, dc, fc);
}

template <class C>
auto for_each_of(C& cpx) {
  return [&cpx](auto&& f) {
    for (auto& fil : cpx.get_filtration_values()) f(fil);
  };
}

BOOST_AUTO_TEST_CASE_TEMPLATE(complex_helpers_bounding_box, Fil, list_of_tested_variants) {
  using T = typename Fil::value_type;
  auto cpx = build_complex<Fil>();

  Box<T> box = get_bounding_box_of_complex(cpx, for_each_of(cpx));
  BOOST_CHECK_EQUAL(box.get_number_of_coordinates(), 2);
  BOOST_CHECK_EQUAL(box.get_lower_corner()[0], 0);
  BOOST_CHECK_EQUAL(box.get_lower_corner()[1], 1);
  BOOST_CHECK_EQUAL(box.get_upper_corner()[0], 4);
  BOOST_CHECK_EQUAL(box.get_upper_corner()[1], 4);

  // empty complex
  Complex_t<Fil> empty;
  Box<T> boxEmpty = get_bounding_box_of_complex(empty, for_each_of(empty));
  BOOST_CHECK(boxEmpty.is_trivial());

  // infinite values are ignored
  auto cpx2 = build_complex<Fil>();
  cpx2.get_filtration_values()[0](0, 0) = Fil::T_inf;
  cpx2.get_filtration_values()[1](0, 1) = Fil::T_m_inf;
  Box<T> box2 = get_bounding_box_of_complex(cpx2, for_each_of(cpx2));
  BOOST_CHECK_EQUAL(box2.get_lower_corner()[0], 1);
  BOOST_CHECK_EQUAL(box2.get_upper_corner()[0], 4);
  BOOST_CHECK_EQUAL(box2.get_lower_corner()[1], 3);
  BOOST_CHECK_EQUAL(box2.get_upper_corner()[1], 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(complex_helpers_normalize, Fil, list_of_tested_variants) {
  using T = typename Fil::value_type;
  using ini = std::initializer_list<T>;

  // automatic bounding box: values are mapped to [0,1]
  auto cpx = build_complex<Fil>();
  normalize_filtration_values_in_complex(cpx, for_each_of(cpx));
  typename Complex_t<Fil>::Filtration_value_container expected = {ini{0, 1}, ini{0.5, 0}, ini{0.25, 2. / 3},
                                                                  ini{0.5, 1}, ini{1, 2. / 3}};
  for (std::size_t i = 0; i < expected.size(); ++i) {
    for (std::size_t p = 0; p < 2; ++p) {
      BOOST_CHECK_CLOSE(cpx.get_filtration_values()[i](0, p) + 1., expected[i](0, p) + 1., 1e-6);
    }
  }

  // given box
  auto cpx2 = build_complex<Fil>();
  Box<T> box({0, 0}, {8, 4});
  normalize_filtration_values_in_complex(cpx2, for_each_of(cpx2), box);
  BOOST_CHECK_CLOSE(cpx2.get_filtration_values()[0](0, 0) + 1., 1., 1e-6);
  BOOST_CHECK_CLOSE(cpx2.get_filtration_values()[0](0, 1), 1., 1e-6);
  BOOST_CHECK_CLOSE(cpx2.get_filtration_values()[4](0, 0), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(cpx2.get_filtration_values()[4](0, 1), 0.75, 1e-6);

  // empty complex: nothing happens
  Complex_t<Fil> empty;
  normalize_filtration_values_in_complex(empty, for_each_of(empty));
  BOOST_CHECK(empty.is_empty());

  // degenerate axis (lower == upper): scale falls back to 1, no division by zero
  Complex_t<Fil> c4(typename Complex_t<Fil>::Boundary_container{{}, {}},
                    typename Complex_t<Fil>::Dimension_container{0, 0},
                    typename Complex_t<Fil>::Filtration_value_container{ini{3, 1}, ini{3, 2}});
  normalize_filtration_values_in_complex(c4, for_each_of(c4));
  BOOST_CHECK_EQUAL(c4.get_filtration_values()[0](0, 0), 0);
  BOOST_CHECK_EQUAL(c4.get_filtration_values()[1](0, 0), 0);
  BOOST_CHECK_EQUAL(c4.get_filtration_values()[0](0, 1), 0);
  BOOST_CHECK_EQUAL(c4.get_filtration_values()[1](0, 1), 1);
}