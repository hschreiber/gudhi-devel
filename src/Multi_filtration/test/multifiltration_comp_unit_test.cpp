/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstdint>      //std::int32_t, std::uint32_t, std::int64_t, std::uint64_t
#include <type_traits>  //std::is_integral_v
#include <utility>      //std::pair
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_filtration"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mp11.hpp>

#include <gudhi/Degree_rips_bifiltration.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>

using Gudhi::multi_filtration::are_equal_filtration_values;
using Gudhi::multi_filtration::Degree_rips_bifiltration;
using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;

struct DR_tag {
  template <typename T>
  using type = Degree_rips_bifiltration<T>;
};

struct DF_tag {
  template <typename T>
  using type = Dynamic_multi_parameter_filtration<T>;
};

struct MF_tag {
  template <typename T>
  using type = Multi_parameter_filtration<T>;
};

using Templates = boost::mp11::mp_list<DR_tag, DF_tag, MF_tag>;
using Types = boost::mp11::mp_list<std::int32_t, std::uint32_t, std::int64_t, std::uint64_t, float, double>;

template <typename Tag, typename T>
using apply_tag = typename Tag::template type<T>;

using FilTypes = boost::mp11::mp_product<apply_tag, Templates, Types>;
using FilPairs = boost::mp11::mp_product<std::pair, FilTypes, FilTypes>;

BOOST_AUTO_TEST_CASE_TEMPLATE(general_fil_equality, FilPair, FilPairs)
{
  using F1 = typename FilPair::first_type;
  using F2 = typename FilPair::second_type;
  using F1_T = typename F1::value_type;
  using F2_T = typename F2::value_type;

  const int numParam = 2;
  std::vector<F1_T> v1 = {5, 0, 6, 1, 3, 2, 4, 3};
  std::vector<F2_T> v2 = {5, 0, 6, 1, 3, 2, 4, 3};

  F1 f1(v1.begin(), v1.end(), numParam);
  F2 f2(v2.begin(), v2.end(), numParam);

  BOOST_CHECK(are_equal_filtration_values(f1, f2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(general_fil_non_equality, FilPair, FilPairs)
{
  using F1 = typename FilPair::first_type;
  using F2 = typename FilPair::second_type;
  using F1_T = typename F1::value_type;
  using F2_T = typename F2::value_type;

  const int numParam = 2;
  std::vector<F1_T> v1 = {5, 0, 6, 1, 3, 2, 4, 3};
  std::vector<F2_T> v2;

  if constexpr (std::is_integral_v<F2_T>) {
    v2 = {5, 0, 6, 1, 4, 2, 3, 3};
  } else {
    v2 = {5.1, 0, 6, 1, 3.6, 2, 4, 3};
  }

  F1 f1(v1.begin(), v1.end(), numParam);
  F2 f2(v2.begin(), v2.end(), numParam);

  BOOST_CHECK(!are_equal_filtration_values(f1, f2));
}
