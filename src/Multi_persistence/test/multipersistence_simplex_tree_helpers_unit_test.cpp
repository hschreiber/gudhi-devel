/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_filtration/Flat_array_filtration.h>
#include <gudhi/Multi_filtration/Degree_bifiltration.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/multi_simplex_tree_helpers.h>

using Gudhi::Simplex_tree;
using Gudhi::multi_filtration::Degree_bifiltration;
using Gudhi::multi_filtration::Flat_array_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration_value;
using Gudhi::multi_persistence::build_degree_filtration;
using Gudhi::multi_persistence::build_simplex_tree_from_complex;
using Gudhi::multi_persistence::fill_axis_with_distance_matrix;
using Gudhi::multi_persistence::fill_axis_with_lowerstar;
using Gudhi::multi_persistence::make_multi_dimensional;
using Gudhi::multi_persistence::make_one_dimensional;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration;

template <typename T>
using Filtration = Multi_parameter_filtration_value<Flat_array_filtration<T>>;
template <typename T>
using Degree_filtration = Multi_parameter_filtration_value<Degree_bifiltration<T>>;
template <typename T>
using Multi_options = Simplex_tree_options_multidimensional_filtration<Filtration<T>>;
template <typename T>
using Std_options = Simplex_tree_options_multidimensional_filtration<T>;
template <typename T>
using Multi_tree = Gudhi::Simplex_tree<Multi_options<T>>;
template <typename T>
using Std_tree = Gudhi::Simplex_tree<Std_options<T>>;
template <typename T>
using Degree_tree = Simplex_tree<Simplex_tree_options_multidimensional_filtration<Degree_filtration<T>>>;
template <typename T>
using Complex = Multi_parameter_filtered_complex<Filtration<T>, std::uint32_t, int>;
template <typename T>
using ini = std::initializer_list<T>;

// "sorted list of vertices" -> filtration value
template <class Tree>
std::map<std::vector<int>, typename Tree::Filtration_value> to_map(Tree& st) {
  std::map<std::vector<int>, typename Tree::Filtration_value> res;
  for (auto sh : st.complex_simplex_range()) {
    std::vector<int> key;
    for (auto v : st.simplex_vertex_range(sh)) key.push_back(v);
    std::sort(key.begin(), key.end());
    res.emplace(key, st.filtration(sh));
  }
  return res;
}

// ---------------------------------------------------------------- make_one_dimensional, make_multi_dimensional

template <typename T>
Multi_tree<T> build_multi_simplex_tree() {
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, ini<T>{0, 1, 2});
  st.insert_simplex_and_subfaces({0, 1, 3}, ini<T>{3, 4, 5});
  st.insert_simplex_and_subfaces({2, 3, 4}, ini<T>{6, 7, 8});
  st.set_num_parameters(3);
  return st;
}

template <typename T>
Std_tree<T> build_std_simplex_tree() {
  Std_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, 1);
  st.insert_simplex_and_subfaces({0, 1, 3}, 4);
  st.insert_simplex_and_subfaces({2, 3, 4}, 7);
  return st;
}

template <typename T>
Multi_tree<T> build_proj_multi_simplex_tree() {
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, ini<T>{0, 1, 2});
  st.insert_simplex_and_subfaces({0, 1, 3}, ini<T>{0, 4, 2});
  st.insert_simplex_and_subfaces({2, 3, 4}, ini<T>{0, 7, 2});
  st.set_num_parameters(3);
  return st;
}

BOOST_AUTO_TEST_CASE(multi_simplex_tree_helpers) {
  Multi_tree<double> multi_st = build_multi_simplex_tree<double>();
  Std_tree<int> std_st = build_std_simplex_tree<int>();
  Multi_tree<float> proj_multi_st = build_proj_multi_simplex_tree<float>();

  Std_tree<int> test1 = make_one_dimensional<Std_options<int>>(multi_st, 1);
  BOOST_CHECK(test1 == std_st);

  Multi_tree<float> test2 = make_multi_dimensional<Multi_options<float>>(std_st, ini<float>{0, 3, 2}, 1);
  BOOST_CHECK(test2 == proj_multi_st);
}

// ---------------------------------------------------------------- build_simplex_tree_from_complex

using T = double;

BOOST_AUTO_TEST_CASE(build_simplex_tree_from_complex_ordered) {
  // 2 vertices and an edge
  Complex<T> cpx({{}, {}, {0, 1}}, {0, 0, 1}, {ini<T>{0, 1}, ini<T>{1, 0}, ini<T>{1, 1}});

  Multi_tree<T> st = build_simplex_tree_from_complex<Multi_options<T>>(cpx);
  BOOST_CHECK_EQUAL(st.num_vertices(), 2);
  BOOST_CHECK_EQUAL(st.num_simplices(), 3);
  BOOST_CHECK_EQUAL(st.dimension(), 1);
  auto m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({0, 1}));
  BOOST_CHECK(m.at({1}) == Filtration<T>({1, 0}));
  BOOST_CHECK(m.at({0, 1}) == Filtration<T>({1, 1}));

  // maxDim
  Multi_tree<T> st0 = build_simplex_tree_from_complex<Multi_options<T>>(cpx, 0);
  BOOST_CHECK_EQUAL(st0.num_simplices(), 2);
  BOOST_CHECK_EQUAL(st0.dimension(), 0);

  // numParam bigger than the current number of parameters: generators are completed with +inf
  Multi_tree<T> st3 = build_simplex_tree_from_complex<Multi_options<T>>(cpx, -1, 3);
  auto m3 = to_map(st3);
  BOOST_CHECK_EQUAL(m3.at({0, 1}).num_parameters(), 3);
  BOOST_CHECK_EQUAL(m3.at({0, 1})(0, 0), 1);
  BOOST_CHECK_EQUAL(m3.at({0, 1})(0, 1), 1);
  BOOST_CHECK_EQUAL(m3.at({0, 1})(0, 2), Filtration<T>::T_inf);

  // numParam smaller: generators are shortened from the end
  Multi_tree<T> st1 = build_simplex_tree_from_complex<Multi_options<T>>(cpx, -1, 1);
  auto m1 = to_map(st1);
  BOOST_CHECK_EQUAL(m1.at({0}).num_parameters(), 1);
  BOOST_CHECK_EQUAL(m1.at({0})(0, 0), 0);
  BOOST_CHECK_EQUAL(m1.at({1})(0, 0), 1);

  // empty complex
  Complex<T> empty;
  Multi_tree<T> stEmpty = build_simplex_tree_from_complex<Multi_options<T>>(empty);
  BOOST_CHECK_EQUAL(stEmpty.num_simplices(), 0);
}

BOOST_AUTO_TEST_CASE(build_simplex_tree_from_complex_not_ordered_by_dimension) {
  // v0, v1, e01, v2, e12
  Complex<T> cpx({{}, {}, {0, 1}, {}, {1, 3}}, {0, 0, 1, 0, 1},
                 {ini<T>{0, 0}, ini<T>{1, 0}, ini<T>{2, 2}, ini<T>{0, 1}, ini<T>{1, 1}});
  BOOST_REQUIRE(!cpx.is_ordered_by_dimension());

  Multi_tree<T> st = build_simplex_tree_from_complex<Multi_options<T>>(cpx);
  BOOST_CHECK_EQUAL(st.num_vertices(), 3);
  BOOST_CHECK_EQUAL(st.num_simplices(), 5);
  auto m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({0, 0}));
  BOOST_CHECK(m.at({1}) == Filtration<T>({1, 0}));
  BOOST_CHECK(m.at({2}) == Filtration<T>({0, 1}));
  BOOST_CHECK(m.at({0, 1}) == Filtration<T>({2, 2}));
  BOOST_CHECK(m.at({1, 2}) == Filtration<T>({1, 1}));

  Multi_tree<T> st0 = build_simplex_tree_from_complex<Multi_options<T>>(cpx, 0);
  BOOST_CHECK_EQUAL(st0.num_simplices(), 3);
}

// ---------------------------------------------------------------- fill_axis_with_lowerstar

BOOST_AUTO_TEST_CASE(fill_axis_with_lowerstar_test) {
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, ini<T>{5, 5});
  st.set_num_parameters(2);

  fill_axis_with_lowerstar(st, std::vector<T>{1, 3, 2}, 1);

  auto m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({5, 1}));
  BOOST_CHECK(m.at({1}) == Filtration<T>({5, 3}));
  BOOST_CHECK(m.at({2}) == Filtration<T>({5, 2}));
  BOOST_CHECK(m.at({0, 1}) == Filtration<T>({5, 3}));
  BOOST_CHECK(m.at({0, 2}) == Filtration<T>({5, 2}));
  BOOST_CHECK(m.at({1, 2}) == Filtration<T>({5, 3}));
  BOOST_CHECK(m.at({0, 1, 2}) == Filtration<T>({5, 3}));

  // other axis, other range value type
  fill_axis_with_lowerstar(st, std::vector<int>{7, 0, 0}, 0);
  m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({7, 1}));
  BOOST_CHECK(m.at({1}) == Filtration<T>({0, 3}));
  BOOST_CHECK(m.at({0, 1, 2}) == Filtration<T>({7, 3}));

  // empty tree
  Multi_tree<T> empty;
  empty.set_num_parameters(2);
  fill_axis_with_lowerstar(empty, std::vector<T>{}, 0);
  BOOST_CHECK_EQUAL(empty.num_simplices(), 0);
}

// ---------------------------------------------------------------- fill_axis_with_distance_matrix

BOOST_AUTO_TEST_CASE(fill_axis_with_distance_matrix_test) {
  Multi_tree<T> st;
  st.insert_simplex_and_subfaces({0, 1, 2}, ini<T>{5, 5});
  st.set_num_parameters(2);

  std::vector<std::vector<T>> dm = {{0, 1, 4}, {1, 0, 2}, {4, 2, 0}};
  fill_axis_with_distance_matrix(st, dm, 0.5, 0);

  auto m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({0.5, 5}));
  BOOST_CHECK(m.at({1}) == Filtration<T>({0.5, 5}));
  BOOST_CHECK(m.at({2}) == Filtration<T>({0.5, 5}));
  BOOST_CHECK(m.at({0, 1}) == Filtration<T>({1, 5}));
  BOOST_CHECK(m.at({0, 2}) == Filtration<T>({4, 5}));
  BOOST_CHECK(m.at({1, 2}) == Filtration<T>({2, 5}));
  BOOST_CHECK(m.at({0, 1, 2}) == Filtration<T>({4, 5}));

  // vertex value bigger than some distances acts as a lower bound
  fill_axis_with_distance_matrix(st, dm, 3, 1);
  m = to_map(st);
  BOOST_CHECK(m.at({0}) == Filtration<T>({0.5, 3}));
  BOOST_CHECK(m.at({0, 1}) == Filtration<T>({1, 3}));
  BOOST_CHECK(m.at({1, 2}) == Filtration<T>({2, 3}));
  BOOST_CHECK(m.at({0, 1, 2}) == Filtration<T>({4, 4}));

  // only the vertices present in the tree are read in the matrix
  Multi_tree<T> st2;
  st2.insert_simplex_and_subfaces({0, 2}, ini<T>{9, 9});
  st2.set_num_parameters(2);
  fill_axis_with_distance_matrix(st2, dm, 0, 0);
  auto m2 = to_map(st2);
  BOOST_CHECK(m2.at({0, 2}) == Filtration<T>({4, 9}));
  BOOST_CHECK(m2.at({0}) == Filtration<T>({0, 9}));
}

// ---------------------------------------------------------------- build_degree_filtration

BOOST_AUTO_TEST_CASE(build_degree_filtration_test) {
  // star with center 0 and leaves 1, 2, 3; edge {0,i} appears at i, vertices at 0
  Std_tree<T> st;
  for (int v = 0; v < 4; ++v) st.insert_simplex({v}, 0.);
  for (int v = 1; v < 4; ++v) st.insert_simplex({0, v}, static_cast<T>(v));

  Degree_tree<T> dt = build_degree_filtration<Degree_tree<T>>(st, std::vector<int>{0, 1, 2});
  BOOST_CHECK_EQUAL(dt.num_vertices(), 4);
  BOOST_CHECK_EQUAL(dt.num_simplices(), 7);

  const T inf = Degree_filtration<T>::T_inf;
  auto m = to_map(dt);
  auto check = [&](const std::vector<int>& s, const std::vector<T>& expected) {
    const auto& f = m.at(s);
    BOOST_REQUIRE_EQUAL(f.num_generators(), expected.size());
    for (std::size_t g = 0; g < expected.size(); ++g) BOOST_CHECK_EQUAL(f(g, 0), expected[g]);
  };
  // center: it has incident cells at 0 (itself), 1, 2, 3
  check({0}, {2, 1, 0});
  // leaves: they have incident cells at 0 (itself) and at the value of their edge only
  check({1}, {inf, 1, 0});
  check({2}, {inf, 2, 0});
  check({3}, {inf, 3, 0});
  // edges: max over the two vertices and of the edge value
  check({0, 1}, {inf, 1, 1});
  check({0, 2}, {inf, 2, 2});
  check({0, 3}, {inf, 3, 3});

  // single degree
  Degree_tree<T> dt1 = build_degree_filtration<Degree_tree<T>>(st, std::vector<int>{1});
  auto m1 = to_map(dt1);
  BOOST_CHECK_EQUAL(m1.at({0}).num_generators(), 1);
  BOOST_CHECK_EQUAL(m1.at({0})(0, 0), 1);
  BOOST_CHECK_EQUAL(m1.at({2})(0, 0), 2);
  BOOST_CHECK_EQUAL(m1.at({0, 3})(0, 0), 3);

  // degree bigger than any vertex degree gives +inf
  Degree_tree<T> dtBig = build_degree_filtration<Degree_tree<T>>(st, std::vector<int>{10});
  auto mBig = to_map(dtBig);
  BOOST_CHECK_EQUAL(mBig.at({0})(0, 0), inf);
  BOOST_CHECK_EQUAL(mBig.at({0, 1})(0, 0), inf);
}
