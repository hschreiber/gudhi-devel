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
#include <map>
#include <stdexcept>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/Slicer.h>
#include <gudhi/Multi_persistence/Persistence_interface_cohomology.h>
#include <gudhi/Multi_persistence/Persistence_interface_homology.h>
#include <gudhi/Multi_persistence/Persistence_interface_vineyard.h>
#include <gudhi/Multi_persistence/tensor_helpers.h>
#include <gudhi/multi_simplex_tree_helpers.h>
#include <gudhi/euler_signed_measure.h>
#include <gudhi/hilbert_signed_measure.h>
#include <gudhi/rank_invariant.h>

using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_persistence::compute_hilbert_signed_measure;
using Gudhi::multi_persistence::compute_rank_invariant;
using Gudhi::multi_persistence::compute_sparse_hilbert_signed_measure;
using Gudhi::multi_persistence::compute_sparse_rank_signed_measure;
using Gudhi::multi_persistence::get_euler_signed_measure;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Persistence_interface_cohomology;
using Gudhi::multi_persistence::Persistence_interface_homology;
using Gudhi::multi_persistence::Persistence_interface_vineyard;
using Gudhi::multi_persistence::Signed_measure;
using Gudhi::multi_persistence::Slicer;

using Filtration = Multi_parameter_filtration_value<Flat_array_filtration<double> >;
using Idx = std::uint32_t;
using D = int;
using Complex = Multi_parameter_filtered_complex<Filtration, Idx, D>;
using ini = std::initializer_list<double>;

struct Multi_persistence_r_options
    : Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::INTRUSIVE_SET, true> {
  using Index = Idx;
  using Dimension = D;
  static const bool has_column_pairings = true;
};

struct Multi_persistence_vineyard_ru_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = true;
};

using list_of_slicers =
    boost::mpl::list<Slicer<Filtration, Persistence_interface_cohomology<Filtration> >,
                     Slicer<Filtration, Persistence_interface_homology<Multi_persistence_r_options, Filtration> >,
                     Slicer<Filtration, Persistence_interface_vineyard<Multi_persistence_vineyard_ru_options> > >;

// two vertices (0,1), (1,0) joined by an edge (1,1)
Complex make_path_complex() { return Complex({{}, {}, {0, 1}}, {0, 0, 1}, {ini{0, 1}, ini{1, 0}, ini{1, 1}}); }

// hollow triangle: 3 vertices, 3 edges
Complex make_cycle_complex() {
  return Complex({{}, {}, {}, {0, 1}, {1, 2}, {0, 2}}, {0, 0, 0, 1, 1, 1},
                 {ini{0, 2}, ini{1, 1}, ini{2, 0}, ini{1, 2}, ini{2, 1}, ini{2, 2}});
}

// Boost macros do not like commas in template arguments
// tests if signed measure is empty
#define EMPTY_SM(...) ((__VA_ARGS__).values.empty())

using Measure = std::map<std::vector<Idx>, int>;

// converts a signed measure into a map "coordinates -> total mass" (order independent)
template <class SM>
Measure to_map(const SM& sm) {
  Measure m;
  const std::size_t n = sm.numberOfParameters;
  BOOST_REQUIRE_EQUAL(sm.coordinates.size(), n * sm.values.size());
  for (std::size_t i = 0; i < sm.values.size(); ++i) {
    std::vector<Idx> key(sm.coordinates.begin() + i * n, sm.coordinates.begin() + (i + 1) * n);
    m[key] += sm.values[i];
  }
  for (auto it = m.begin(); it != m.end();) it = (it->second == 0) ? m.erase(it) : std::next(it);
  return m;
}

// ------------------------------------------------------------------ Euler

BOOST_AUTO_TEST_CASE(euler_signed_measure) {
  using Opt = Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration<Filtration>;
  Gudhi::Simplex_tree<Opt> st;
  st.insert_simplex({0}, ini{0, 1});
  st.insert_simplex({1}, ini{1, 0});
  st.insert_simplex({0, 1}, ini{1, 1});
  st.set_num_parameters(2);
  const std::vector<int> shape = {3, 3};

  // Mobius inversion, no padding: signed point masses at the filtration values
  auto sm = get_euler_signed_measure<int, Idx>(st, shape, false);
  BOOST_CHECK_EQUAL(sm.numberOfParameters, 2);
  BOOST_CHECK((to_map(sm) == Measure{{{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, -1}}));
  // default value of mobiusInversion is true
  BOOST_CHECK((to_map(get_euler_signed_measure<int, Idx>(st, shape, false, true)) == to_map(sm)));
  // default template arguments
  auto smDefault = get_euler_signed_measure(st, shape, false);
  BOOST_CHECK_EQUAL(smDefault.values.size(), 3);

  // Mobius inversion with zero padding: the last grid index is only used as padding
  auto smPad = get_euler_signed_measure<int, Idx>(st, shape, true);
  BOOST_CHECK((to_map(smPad) ==
               Measure{{{0, 1}, 1}, {{0, 2}, -1}, {{1, 0}, 1}, {{1, 1}, -1}, {{2, 0}, -1}, {{2, 2}, 1}}));

  // no Mobius inversion: the Euler characteristic surface itself (1 everywhere except where nothing is born)
  auto smSurface = get_euler_signed_measure<int, Idx>(st, shape, false, false);
  Measure expectedSurface;
  for (Idx i = 0; i < 3; ++i) {
    for (Idx j = 0; j < 3; ++j) {
      if (i != 0 || j != 0) expectedSurface[{i, j}] = 1;
    }
  }
  BOOST_CHECK(to_map(smSurface) == expectedSurface);

  // verbose does not change the result
  BOOST_CHECK(to_map(get_euler_signed_measure<int, Idx>(st, shape, false, true, true)) == to_map(sm));

  // empty grid
  auto smEmpty = get_euler_signed_measure<int, Idx>(st, std::vector<int>{0, 3}, false);
  BOOST_CHECK(smEmpty.values.empty());
  BOOST_CHECK(smEmpty.coordinates.empty());

  // empty complex
  Gudhi::Simplex_tree<Opt> stEmpty;
  stEmpty.set_num_parameters(2);
  BOOST_CHECK(EMPTY_SM(get_euler_signed_measure<int, Idx>(stEmpty, shape, false)));

  // a simplex of dimension 2: 3 vertices - 3 edges + 1 triangle = 1, all born at (0,0)
  Gudhi::Simplex_tree<Opt> st2;
  st2.insert_simplex_and_subfaces({0, 1, 2}, ini{0, 0});
  st2.set_num_parameters(2);
  BOOST_CHECK((to_map(get_euler_signed_measure<int, Idx>(st2, shape, false)) == Measure{{{0, 0}, 1}}));
}

// ------------------------------------------------------------------ Hilbert

BOOST_AUTO_TEST_CASE_TEMPLATE(hilbert_signed_measure, Slicer_t, list_of_slicers) {
  const std::vector<int> shape = {3, 3};
  const std::vector<int> degrees = {0};

  {
    Slicer_t slicerDense(make_path_complex());
    auto sm = compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, degrees);
    // first coordinate of each point is the index in degrees
    BOOST_CHECK_EQUAL(sm.numberOfParameters, 3);
    BOOST_CHECK((to_map(sm) == Measure{{{0, 0, 1}, 1}, {{0, 1, 0}, 1}, {{0, 1, 1}, -1}}));

    Slicer_t slicerSparse(make_path_complex());
    auto sms = compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, degrees);
    BOOST_CHECK_EQUAL(sms.numberOfParameters, 3);
    BOOST_CHECK(to_map(sms) == to_map(sm));
  }
  {
    // zero padding
    Slicer_t slicerDense(make_path_complex());
    auto sm = compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, degrees, true);
    Measure expected = {{{0, 0, 1}, 1}, {{0, 0, 2}, -1}, {{0, 1, 1}, -1}, {{0, 1, 2}, 1}};
    BOOST_CHECK(to_map(sm) == expected);
    Slicer_t slicerSparse(make_path_complex());
    BOOST_CHECK(to_map(compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, degrees, true)) ==
                expected);
  }
  {
    // all optional arguments
    Slicer_t slicerDense(make_path_complex());
    auto sm = compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, degrees, false, false, true, 1);
    BOOST_CHECK(to_map(sm) == (Measure{{{0, 0, 1}, 1}, {{0, 1, 0}, 1}, {{0, 1, 1}, -1}}));
    Slicer_t slicerSparse(make_path_complex());
    auto sms = compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, degrees, false, false, true, 1);
    BOOST_CHECK(to_map(sms) == to_map(sm));
  }
  {
    // degree 1 of a tree is trivial
    Slicer_t slicerDense(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, std::vector<int>{1})));
    Slicer_t slicerSparse(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, std::vector<int>{1})));
  }
  {
    // empty degrees
    Slicer_t slicerDense(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, std::vector<int>{})));
    Slicer_t slicerSparse(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, std::vector<int>{})));
  }
  {
    // empty grid
    Slicer_t slicerDense(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_hilbert_signed_measure<int, Idx>(slicerDense, std::vector<int>{0, 3}, degrees)));
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(hilbert_signed_measure_cycle, Slicer_t, list_of_slicers) {
  // dense and sparse versions must agree, for several degrees at once
  const std::vector<int> shape = {4, 4};  // the last index is only padding when zeroPad is true
  const std::vector<int> degrees = {0, 1};
  for (bool zeroPad : {false, true}) {
    Slicer_t slicerDense(make_cycle_complex());
    auto dense = compute_hilbert_signed_measure<int, Idx>(slicerDense, shape, degrees, zeroPad);
    Slicer_t slicerSparse(make_cycle_complex());
    auto sparse = compute_sparse_hilbert_signed_measure<int, Idx>(slicerSparse, shape, degrees, zeroPad);
    BOOST_CHECK(!dense.values.empty());
    BOOST_CHECK(to_map(dense) == to_map(sparse));
  }

  // without padding: the cycle is born at (2,2)
  Slicer_t slicer(make_cycle_complex());
  auto sm = compute_hilbert_signed_measure<int, Idx>(slicer, shape, std::vector<int>{1});
  BOOST_CHECK((to_map(sm) == Measure{{{0, 2, 2}, 1}}));
  // H0: components born at (0,2), (1,1) and (2,0), merged pairwise at (1,2) and (2,1)
  Slicer_t slicerH0(make_cycle_complex());
  auto smH0 = compute_hilbert_signed_measure<int, Idx>(slicerH0, shape, std::vector<int>{0});
  BOOST_CHECK((to_map(smH0) ==
               Measure{{{0, 0, 2}, 1}, {{0, 1, 1}, 1}, {{0, 2, 0}, 1}, {{0, 1, 2}, -1}, {{0, 2, 1}, -1}}));
}

BOOST_AUTO_TEST_CASE(hilbert_signed_measure_unsupported_number_of_parameters) {
  using Slicer_t = Slicer<Filtration, Persistence_interface_cohomology<Filtration> >;
  Complex cpx({{}, {}}, {0, 0}, {ini{0, 0, 0, 0, 0}, ini{1, 1, 1, 1, 1}});
  Slicer_t slicer(cpx);
  BOOST_CHECK_THROW(
      (compute_sparse_hilbert_signed_measure<int, Idx>(slicer, std::vector<int>{2, 2, 2, 2, 2}, std::vector<int>{0})),
      std::runtime_error);
}

// ------------------------------------------------------------------ Rank invariant

BOOST_AUTO_TEST_CASE_TEMPLATE(sparse_rank_signed_measure, Slicer_t, list_of_slicers) {
  const std::vector<Idx> shape = {3, 3};
  const std::vector<Idx> degrees = {0};

  // for is_vine, the ignoreInf parameter is ignored so it is useless to test it at true
  if constexpr (!Slicer_t::Persistence::is_vine) {
    Slicer_t slicer(make_path_complex());
    auto sm = compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, degrees, false, true);
    // (degree index, birth, death) with death = shape meaning "never"
    BOOST_CHECK_EQUAL(sm.numberOfParameters, 5);
    BOOST_CHECK((to_map(sm) == Measure{{{0, 0, 1, 3, 3}, 1}, {{0, 1, 0, 3, 3}, 1}, {{0, 1, 1, 3, 3}, -1}}));
  }
  if constexpr (!Slicer_t::Persistence::is_vine) {
    Slicer_t slicer(make_path_complex());
    auto sm = compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, degrees, true, true);
    BOOST_CHECK((to_map(sm) == Measure{{{0, 0, 1, 2, 2}, 1}, {{0, 1, 0, 2, 2}, 1}, {{0, 1, 1, 2, 2}, -1}}));
  }
  {
    // ignoreInf = false, n_jobs = 1
    Slicer_t slicer(make_path_complex());
    auto sm = compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, degrees, false, false, 1);
    BOOST_CHECK_EQUAL(sm.numberOfParameters, 5);
    BOOST_CHECK(!sm.values.empty());
  }
  {
    Slicer_t slicer(make_path_complex());
    BOOST_CHECK(EMPTY_SM(compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, std::vector<Idx>{}, false, true)));
    BOOST_CHECK(EMPTY_SM(compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, std::vector<Idx>{1}, false, true)));
  }
  if constexpr (!Slicer_t::Persistence::is_vine) {
    // H1 of the cycle: born at (2,2), never dies
    Slicer_t slicer(make_cycle_complex());
    auto sm = compute_sparse_rank_signed_measure<int, Idx>(slicer, shape, std::vector<Idx>{1}, false, true);
    BOOST_CHECK((to_map(sm) == Measure{{{0, 2, 2, 3, 3}, 1}}));
  }
  {
    // three parameters
    Complex cpx({{}, {}, {0, 1}}, {0, 0, 1}, {ini{0, 1, 0}, ini{1, 0, 0}, ini{1, 1, 1}});
    Slicer_t slicer(cpx);
    auto sm = compute_sparse_rank_signed_measure<int, Idx>(slicer, std::vector<Idx>{2, 2, 2}, degrees, false, true);
    BOOST_CHECK_EQUAL(sm.numberOfParameters, 7);
    BOOST_CHECK(!sm.values.empty());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(dense_rank_invariant, Slicer_t, list_of_slicers) {
  const std::vector<Idx> shape = {3, 3};
  const std::vector<Idx> degrees = {0};

  Slicer_t slicer(make_path_complex());
  auto out = compute_rank_invariant<int, Idx>(slicer, shape, degrees, false);
  const std::vector<int>& values = out.first;
  const auto& tensorShape = out.second;

  // shape is (degrees, birth coordinates, death coordinates)
  BOOST_REQUIRE_EQUAL(tensorShape.size(), 5);
  BOOST_CHECK_EQUAL(tensorShape[0], 1);
  for (std::size_t i = 1; i < 5; ++i) BOOST_CHECK_EQUAL(tensorShape[i], 3);
  BOOST_REQUIRE_EQUAL(values.size(), 81);

  // rank(H_0(a) -> H_0(b)) is 1 iff a <= b and a != (0,0), and 0 otherwise
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          const bool expected = (i <= k && j <= l && (i != 0 || j != 0));
          BOOST_CHECK_EQUAL(values[((i * 3 + j) * 3 + k) * 3 + l], expected ? 1 : 0);
        }
      }
    }
  }

  // empty degrees: nothing to compute
  Slicer_t slicer2(make_path_complex());
  auto outEmpty = compute_rank_invariant<int, Idx>(slicer2, shape, std::vector<Idx>{}, false);
  BOOST_CHECK(outEmpty.first.empty());
}

// The sparse signed measure integrates to the dense rank invariant:
// rank(a -> b) = sum of the masses of the points (a', b') with a' <= a and b <= b'
// (a death coordinate equal to the extent of the grid means "never dies").
template <class Slicer_t>
void check_sparse_rank_is_consistent_with_dense(bool ignoreInf) {
  const std::vector<Idx> shape = {3, 3};
  const std::vector<Idx> degrees = {0, 1};

  Slicer_t slicerDense(make_cycle_complex());
  auto dense = compute_rank_invariant<int, Idx>(slicerDense, shape, degrees, ignoreInf);
  Slicer_t slicerSparse(make_cycle_complex());
  auto sparse = compute_sparse_rank_signed_measure<int, Idx>(slicerSparse, shape, degrees, false, ignoreInf);
  BOOST_REQUIRE_EQUAL(sparse.numberOfParameters, 5);
  BOOST_REQUIRE_EQUAL(dense.first.size(), 2 * 81);

  for (Idx d = 0; d < 2; ++d) {
    for (Idx i = 0; i < 3; ++i) {
      for (Idx j = 0; j < 3; ++j) {
        for (Idx k = i; k < 3; ++k) {
          for (Idx l = j; l < 3; ++l) {
            int integrated = 0;
            for (std::size_t p = 0; p < sparse.values.size(); ++p) {
              const Idx* c = &sparse.coordinates[5 * p];
              if (c[0] == d && c[1] <= i && c[2] <= j && k <= c[3] && l <= c[4]) integrated += sparse.values[p];
            }
            BOOST_CHECK_EQUAL(integrated, dense.first[(((d * 3 + i) * 3 + j) * 3 + k) * 3 + l]);
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sparse_rank_consistent_with_dense_ignore_inf, Slicer_t, list_of_slicers) {
  check_sparse_rank_is_consistent_with_dense<Slicer_t>(true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sparse_rank_consistent_with_dense_keep_inf, Slicer_t, list_of_slicers) {
  check_sparse_rank_is_consistent_with_dense<Slicer_t>(false);
}

