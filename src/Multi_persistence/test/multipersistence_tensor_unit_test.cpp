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
#include <limits>
#include <stdexcept>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/tensor_helpers.h>

using Gudhi::multi_persistence::differentiate;
using Gudhi::multi_persistence::Signed_measure;
using Gudhi::multi_persistence::sparsify;
using Gudhi::multi_persistence::Tensor;

BOOST_AUTO_TEST_CASE(tensor_constructors_and_access) {
  std::vector<int> data(24, 0);
  Tensor<int> t(data, std::vector<std::size_t>{2, 3, 4});

  BOOST_CHECK_EQUAL(t.rank(), 3);
  BOOST_CHECK_EQUAL(t.rank_dynamic(), 3);
  BOOST_CHECK_EQUAL(t.size(), 24);
  BOOST_CHECK(!t.empty());
  BOOST_CHECK_EQUAL(t.extent(0), 2);
  BOOST_CHECK_EQUAL(t.extent(1), 3);
  BOOST_CHECK_EQUAL(t.extent(2), 4);
  BOOST_CHECK_EQUAL(t.stride(0), 12);
  BOOST_CHECK_EQUAL(t.stride(1), 4);
  BOOST_CHECK_EQUAL(t.stride(2), 1);
  BOOST_CHECK(t.extents() == std::vector<std::size_t>({2, 3, 4}));
  BOOST_CHECK_EQUAL(t.data_handle(), data.data());
  BOOST_CHECK_EQUAL(Tensor<int>::static_extent(0), std::numeric_limits<std::size_t>::max());

  BOOST_CHECK(Tensor<int>::is_always_unique());
  BOOST_CHECK(Tensor<int>::is_always_exhaustive());
  BOOST_CHECK(Tensor<int>::is_always_strided());
  BOOST_CHECK(t.is_unique());
  BOOST_CHECK(t.is_exhaustive());
  BOOST_CHECK(t.is_strided());

  t(1, 2, 3) = 7;
  BOOST_CHECK_EQUAL(data[23], 7);
  t(0, 1, 2) = 5;
  BOOST_CHECK_EQUAL(data[6], 5);
  BOOST_CHECK_EQUAL((t[std::vector<std::size_t>{1, 2, 3}]), 7);
  BOOST_CHECK_EQUAL((t[{0, 1, 2}]), 5);
  BOOST_CHECK_EQUAL(t(0, 0, 0), 0);

  std::vector<int> data2(6, 1);
  Tensor<int> t2(data2, std::vector<std::size_t>{3, 2});
  BOOST_CHECK_EQUAL(t2.rank(), 2);
  BOOST_CHECK_EQUAL(t2.size(), 6);

  // initializer-list shape constructor
  Tensor<int> t3(data2, std::initializer_list<int>{2, 3});
  BOOST_CHECK_EQUAL(t3.extent(1), 3);

  BOOST_CHECK_THROW(t.extent(3), std::out_of_range);
  BOOST_CHECK_THROW(t.stride(3), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(tensor_swap) {
  std::vector<int> d1(6, 1), d2(4, 2);
  Tensor<int> a(d1, std::vector<std::size_t>{2, 3});
  Tensor<int> b(d2, std::vector<std::size_t>{4});
  swap(a, b);
  BOOST_CHECK_EQUAL(a.rank(), 1);
  BOOST_CHECK_EQUAL(a.size(), 4);
  BOOST_CHECK_EQUAL(a.data_handle(), d2.data());
  BOOST_CHECK_EQUAL(b.rank(), 2);
  BOOST_CHECK_EQUAL(b.size(), 6);
  BOOST_CHECK_EQUAL(b.data_handle(), d1.data());
}

BOOST_AUTO_TEST_CASE(tensor_differentiate) {
  // 1D
  std::vector<int> d1 = {1, 3, 6, 10};
  Tensor<int> t1(d1, std::vector<std::size_t>{4});
  differentiate(t1, 0);
  BOOST_CHECK((d1 == std::vector<int>{1, 2, 3, 4}));

  // 2D, shape (2,3), row-major
  std::vector<int> d2 = {1, 2, 3, 5, 7, 9};
  Tensor<int> t2(d2, std::vector<std::size_t>{2, 3});
  differentiate(t2, 0);
  BOOST_CHECK((d2 == std::vector<int>{1, 2, 3, 4, 5, 6}));

  std::vector<int> d3 = {1, 3, 6, 4, 9, 15};
  Tensor<int> t3(d3, std::vector<std::size_t>{2, 3});
  differentiate(t3, 1);
  BOOST_CHECK((d3 == std::vector<int>{1, 2, 3, 4, 5, 6}));

  // 3D, axis in the middle
  std::vector<int> d4 = {1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3};
  Tensor<int> t4(d4, std::vector<std::size_t>{2, 3, 2});
  differentiate(t4, 1);
  BOOST_CHECK((d4 == std::vector<int>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

  // axis with extent 1: no-op
  std::vector<int> d5 = {4, 5, 6};
  Tensor<int> t5(d5, std::vector<std::size_t>{1, 3});
  differentiate(t5, 0);
  BOOST_CHECK((d5 == std::vector<int>{4, 5, 6}));

  // n_jobs argument
  std::vector<int> d6 = {1, 3, 6, 10};
  Tensor<int> t6(d6, std::vector<std::size_t>{4});
  differentiate(t6, 0, 1);
  BOOST_CHECK((d6 == std::vector<int>{1, 2, 3, 4}));

  // large tensor (goes through the parallel branch when TBB is used)
  const std::size_t n = 1u << 16;
  std::vector<long> big(n);
  for (std::size_t i = 0; i < n; ++i) big[i] = static_cast<long>(i * (i + 1) / 2);
  Tensor<long> tb(big, std::vector<std::size_t>{256, 256});
  differentiate(tb, 1);
  for (std::size_t r = 0; r < 256; ++r) {
    for (std::size_t c = 1; c < 256; ++c) {
      // with require instead of check, the loop stops at the first error avoiding the print of 256 * 256 errors...
      BOOST_REQUIRE_EQUAL(big[r * 256 + c], static_cast<long>(r * 256 + c));
    }
  }
}

BOOST_AUTO_TEST_CASE(tensor_sparsify) {
  std::vector<int> d = {0, 2, 0, 0, 0, -3};
  Tensor<int> t(d, std::vector<std::size_t>{2, 3});

  Signed_measure<std::uint32_t, int> sm = sparsify<std::uint32_t>(t);
  BOOST_CHECK_EQUAL(sm.numberOfParameters, 2);
  BOOST_CHECK((sm.coordinates == std::vector<std::uint32_t>{0, 1, 1, 2}));
  BOOST_CHECK((sm.values == std::vector<int>{2, -3}));

  // flip the second axis: coordinate becomes extent - value
  auto sm2 = sparsify<std::uint32_t>(t, {false, true});
  BOOST_CHECK_EQUAL(sm2.numberOfParameters, 2);
  BOOST_CHECK((sm2.coordinates == std::vector<std::uint32_t>{0, 2, 1, 1}));
  BOOST_CHECK((sm2.values == std::vector<int>{2, -3}));

  // flip first axis only
  auto sm3 = sparsify<std::uint32_t>(t, {true});
  BOOST_CHECK((sm3.coordinates == std::vector<std::uint32_t>{2, 1, 1, 2}));

  // zero tensor
  std::vector<int> z(6, 0);
  Tensor<int> tz(z, std::vector<std::size_t>{3, 2});
  auto sm4 = sparsify<std::uint16_t>(tz);
  BOOST_CHECK_EQUAL(sm4.numberOfParameters, 2);
  BOOST_CHECK(sm4.coordinates.empty());
  BOOST_CHECK(sm4.values.empty());
}
