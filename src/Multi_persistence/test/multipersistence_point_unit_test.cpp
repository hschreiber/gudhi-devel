/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Point.h>

using Gudhi::multi_persistence::Point;

using list_of_tested_variants = boost::mpl::list<double, float, int>;

BOOST_AUTO_TEST_CASE_TEMPLATE(point_constructors, T, list_of_tested_variants) {
  Point<T> p0;
  BOOST_CHECK_EQUAL(p0.size(), 0);

  Point<T> p1(3);
  BOOST_CHECK_EQUAL(p1.size(), 3);

  Point<T> p2(3, 5);
  BOOST_CHECK_EQUAL(p2.size(), 3);
  for (std::size_t i = 0; i < 3; ++i) BOOST_CHECK_EQUAL(p2[i], 5);

  std::vector<T> v = {1, 2, 3, 4};
  Point<T> p3(v.begin(), v.end());
  BOOST_CHECK_EQUAL(p3.size(), 4);
  for (std::size_t i = 0; i < 4; ++i) BOOST_CHECK_EQUAL(p3[i], v[i]);

  Point<T> p4 = {1, 2, 3, 4};
  BOOST_CHECK(p3 == p4);

  Point<T> p5(v);  // copy of container
  BOOST_CHECK(p5 == p4);

  std::vector<T> v2 = v;
  Point<T> p6(std::move(v2));  // move of container
  BOOST_CHECK(p6 == p4);

  Point<T> p7(p4);  // copy
  BOOST_CHECK(p7 == p4);

  Point<T> p8(std::move(p7));  // move
  BOOST_CHECK(p8 == p4);

  Point<T> p9;
  p9 = {7, 8};  // assignment from initializer list
  BOOST_CHECK_EQUAL(p9.size(), 2);
  BOOST_CHECK_EQUAL(p9[0], 7);
  BOOST_CHECK_EQUAL(p9[1], 8);

  Point<T> p10;
  p10 = p4;  // copy assignment
  BOOST_CHECK(p10 == p4);
  Point<T> p11;
  p11 = std::move(p10);  // move assignment
  BOOST_CHECK(p11 == p4);

  std::vector<T> converted = p4;  // conversion operator
  BOOST_CHECK(converted == v);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_access, T, list_of_tested_variants) {
  Point<T> p = {1, 2, 3};
  const Point<T>& cp = p;

  BOOST_CHECK_EQUAL(p.size(), 3);
  BOOST_CHECK_EQUAL(p[0], 1);
  BOOST_CHECK_EQUAL(cp[1], 2);
  BOOST_CHECK_EQUAL(p.at(2), 3);
  BOOST_CHECK_EQUAL(cp.at(2), 3);
  BOOST_CHECK_THROW(p.at(3), std::out_of_range);
  BOOST_CHECK_THROW(cp.at(3), std::out_of_range);
  BOOST_CHECK_EQUAL(p.front(), 1);
  BOOST_CHECK_EQUAL(cp.front(), 1);
  BOOST_CHECK_EQUAL(p.back(), 3);
  BOOST_CHECK_EQUAL(cp.back(), 3);

  p[0] = 10;
  p.at(1) = 20;
  p.front() += 1;
  p.back() = 30;
  BOOST_CHECK_EQUAL(p[0], 11);
  BOOST_CHECK_EQUAL(p[1], 20);
  BOOST_CHECK_EQUAL(p[2], 30);

  std::vector<T>& container = p.retrieve_underlying_container();
  BOOST_CHECK_EQUAL(container.size(), 3);
  container.push_back(40);
  BOOST_CHECK_EQUAL(p.size(), 4);
  BOOST_CHECK_EQUAL(cp.retrieve_underlying_container()[3], 40);

  BOOST_CHECK_EQUAL(p.data(), container.data());
  BOOST_CHECK_EQUAL(cp.data(), container.data());

  p.resize(2);
  BOOST_CHECK_EQUAL(p.size(), 2);
  p.resize(0);
  BOOST_CHECK_EQUAL(p.size(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_iterators, T, list_of_tested_variants) {
  Point<T> p = {1, 2, 3};
  const Point<T>& cp = p;

  BOOST_CHECK_EQUAL(std::distance(p.begin(), p.end()), 3);
  BOOST_CHECK_EQUAL(std::distance(cp.begin(), cp.end()), 3);
  BOOST_CHECK_EQUAL(std::distance(p.cbegin(), p.cend()), 3);

  BOOST_CHECK_EQUAL(*p.begin(), 1);
  BOOST_CHECK_EQUAL(*cp.begin(), 1);
  BOOST_CHECK_EQUAL(*p.cbegin(), 1);
  BOOST_CHECK_EQUAL(*(p.end() - 1), 3);
  BOOST_CHECK_EQUAL(*(cp.end() - 1), 3);
  BOOST_CHECK_EQUAL(*(p.cend() - 1), 3);

  BOOST_CHECK_EQUAL(*p.rbegin(), 3);
  BOOST_CHECK_EQUAL(*cp.rbegin(), 3);
  BOOST_CHECK_EQUAL(*p.crbegin(), 3);
  BOOST_CHECK_EQUAL(*(p.rend() - 1), 1);
  BOOST_CHECK_EQUAL(*(cp.rend() - 1), 1);
  BOOST_CHECK_EQUAL(*(p.crend() - 1), 1);

  for (auto it = p.begin(); it != p.end(); ++it) *it += 1;
  BOOST_CHECK_EQUAL(p[0], 2);
  BOOST_CHECK_EQUAL(p[1], 3);
  BOOST_CHECK_EQUAL(p[2], 4);

  T sum = 0;
  for (T v : cp) sum += v;
  BOOST_CHECK_EQUAL(sum, 9);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_swap, T, list_of_tested_variants) {
  Point<T> a = {1, 2, 3};
  Point<T> b = {4, 5};

  a.swap(b);
  BOOST_CHECK(a == Point<T>({4, 5}));
  BOOST_CHECK(b == Point<T>({1, 2, 3}));

  swap(a, b);
  BOOST_CHECK(a == Point<T>({1, 2, 3}));
  BOOST_CHECK(b == Point<T>({4, 5}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_comparisons, T, list_of_tested_variants) {
  Point<T> a = {1, 2, 3};
  Point<T> b = {1, 2, 3};
  Point<T> c = {2, 3, 4};
  Point<T> d = {2, 1, 4};  // incomparable with a and c
  Point<T> e = {1, 2};

  BOOST_CHECK(a == a);
  BOOST_CHECK(a == b);
  BOOST_CHECK(!(a != b));
  BOOST_CHECK(a != c);
  BOOST_CHECK(a != e);  // different sizes
  BOOST_CHECK(!(a == e));

  BOOST_CHECK(a <= a);
  BOOST_CHECK(a <= b);
  BOOST_CHECK(a <= c);
  BOOST_CHECK(!(c <= a));
  BOOST_CHECK(!(a <= d));
  BOOST_CHECK(!(d <= a));

  BOOST_CHECK(!(a < a));
  BOOST_CHECK(!(a < b));
  BOOST_CHECK(a < c);
  BOOST_CHECK(!(c < a));
  BOOST_CHECK(!(a < d));
  BOOST_CHECK(!(d < a));

  BOOST_CHECK(c >= a);
  BOOST_CHECK(a >= b);
  BOOST_CHECK(!(a >= c));
  BOOST_CHECK(c > a);
  BOOST_CHECK(!(a > b));
  BOOST_CHECK(!(a > c));
  BOOST_CHECK(!(d >= a));
  BOOST_CHECK(!(d > a));

  Point<T> inf(3, Point<T>::T_inf);
  Point<T> minf(3, Point<T>::T_m_inf);
  BOOST_CHECK(minf <= a);
  BOOST_CHECK(a <= inf);
  BOOST_CHECK(minf < a);
  BOOST_CHECK(a < inf);
}

BOOST_AUTO_TEST_CASE(point_comparisons_nan) {
  using P = Point<double>;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  P a = {1, nan, 3};
  P b = {1, 2, 3};
  BOOST_CHECK(!(a < b));
  BOOST_CHECK(!(b < a));
  BOOST_CHECK(!(a <= b));
  BOOST_CHECK(!(b <= a));
  BOOST_CHECK(!(a >= b));
  BOOST_CHECK(!(a > b));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_arithmetic_with_scalar, T, list_of_tested_variants) {
  Point<T> p = {2, 4, 6};

  BOOST_CHECK(-p == Point<T>({-2, -4, -6}));

  BOOST_CHECK(p + T(1) == Point<T>({3, 5, 7}));
  BOOST_CHECK(T(1) + p == Point<T>({3, 5, 7}));
  BOOST_CHECK(p - T(1) == Point<T>({1, 3, 5}));
  BOOST_CHECK(T(10) - p == Point<T>({8, 6, 4}));
  BOOST_CHECK(p * T(2) == Point<T>({4, 8, 12}));
  BOOST_CHECK(T(2) * p == Point<T>({4, 8, 12}));
  BOOST_CHECK(p / T(2) == Point<T>({1, 2, 3}));
  BOOST_CHECK(T(12) / p == Point<T>({6, 3, 2}));

  Point<T> q = p;
  (q += T(1)) += T(1);
  BOOST_CHECK(q == Point<T>({4, 6, 8}));
  (q -= T(1)) -= T(1);
  BOOST_CHECK(q == p);
  (q *= T(2)) *= T(2);
  BOOST_CHECK(q == Point<T>({8, 16, 24}));
  (q /= T(2)) /= T(2);
  BOOST_CHECK(q == p);

  // the binary operators must not modify the original
  BOOST_CHECK(p == Point<T>({2, 4, 6}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_arithmetic_with_point, T, list_of_tested_variants) {
  Point<T> p = {2, 4, 6};
  Point<T> q = {1, 2, 3};

  BOOST_CHECK(p + q == Point<T>({3, 6, 9}));
  BOOST_CHECK(p - q == Point<T>({1, 2, 3}));
  BOOST_CHECK(p * q == Point<T>({2, 8, 18}));
  BOOST_CHECK(p / q == Point<T>({2, 2, 2}));

  Point<T> r = p;
  r += q;
  BOOST_CHECK(r == Point<T>({3, 6, 9}));
  r -= q;
  BOOST_CHECK(r == p);
  r *= q;
  BOOST_CHECK(r == Point<T>({2, 8, 18}));
  r /= q;
  BOOST_CHECK(r == p);

  // operations between points of different value types
  Point<int> pi = {1, 1, 1};
  Point<T> s = {1, 2, 3};
  BOOST_CHECK(s + pi == Point<T>({2, 3, 4}));
  BOOST_CHECK(s - pi == Point<T>({0, 1, 2}));
  BOOST_CHECK(s * pi == s);
  BOOST_CHECK(s / pi == s);
}

BOOST_AUTO_TEST_CASE(point_arithmetic_infinity) {
  using P = Point<double>;
  const double inf = P::T_inf;
  const double minf = P::T_m_inf;

  P p = {inf, minf, 1};
  BOOST_CHECK(p + 1.0 == P({inf, minf, 2}));
  BOOST_CHECK(p - 1.0 == P({inf, minf, 0}));
  BOOST_CHECK(-p == P({minf, inf, -1}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(point_serialization, T, list_of_tested_variants) {
  Point<T> p = {1, 2, 3, 4};
  const std::size_t size = get_serialization_size_of(p);
  BOOST_CHECK(size > 0);

  std::vector<char> buffer(size);
  char* end = serialize_value_to_char_buffer(p, buffer.data());
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(end - buffer.data()), size);

  Point<T> q;
  const char* cend = deserialize_value_from_char_buffer(q, buffer.data());
  BOOST_CHECK_EQUAL(static_cast<std::size_t>(cend - buffer.data()), size);
  BOOST_CHECK(p == q);
}
