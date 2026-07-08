/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file multi_filtration_conversions.h
 * @author Hannah Schreiber
 */

#ifndef MF_CONVERSIONS_H_
#define MF_CONVERSIONS_H_

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <limits>
#include <cmath>

#include <boost/safe_numerics/safe_integer.hpp>

#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Degree_rips_bifiltration.h>

namespace Gudhi {

namespace multi_filtration {

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given multi filtration value into the type given as template argument. It is assumed that the
 * given value is simplified (i.e. minimal and ordered lexicographically). If the new type is
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration it is additionally assumed that the given value is compatible
 *  with the type, that is,the number of parameters is 2 and the second parameter is an index (positive and convertible
 * to an integer without loss).
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
inline Out_multi_filtration as_type(const Multi_parameter_filtration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
    return Out_multi_filtration(f.begin(), f.end(), f.num_parameters());
  } else if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");
    U inf = co ? Out_multi_filtration::T_m_inf : Out_multi_filtration::T_inf;

    T maxIndex = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      maxIndex = maxIndex < f(g, 1) ? f(g, 1) : maxIndex;
    }

    std::vector<U> values(maxIndex + 1, inf);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 1)] = f(g, 0);
    }
    return Out_multi_filtration(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given multi filtration value into the type given as template argument. It is assumed that the
 * given value is simplified (i.e. minimal and ordered lexicographically). If the new type is
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration it is additionally assumed that the given value is compatible
 * with the type, that is,the number of parameters is 2 and the second parameter is an index (positive and convertible
 * to an integer without loss).
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
inline Out_multi_filtration as_type(const Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
    std::vector<U> values(f.num_entries());
    std::size_t i = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      for (std::size_t p = 0; p < f.num_parameters(); ++p) {
        values[i] = f(g, p);
        ++i;
      }
    }
    return Out_multi_filtration(std::move(values), f.num_parameters());
  } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    if (f.num_parameters() != 2) throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");
    U inf = co ? Out_multi_filtration::T_m_inf : Out_multi_filtration::T_inf;

    T maxIndex = 0;
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      maxIndex = maxIndex < f(g, 1) ? f(g, 1) : maxIndex;
    }

    std::vector<U> values(maxIndex + 1, inf);
    for (std::size_t g = 0; g < f.num_generators(); ++g) {
      values[f(g, 1)] = f(g, 0);
    }
    return Out_multi_filtration(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given multi filtration value into the type given as template argument. If the new type is
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration it is additionally assumed that the given value is compatible
 * with the type, that is,the number of parameters is 2 and the second parameter is an index (positive and convertible
 * to an integer without loss).
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 * @param numberOfParameters Number of parameters of the filtration value.
 */
template <class Out_multi_filtration, typename T>
inline Out_multi_filtration as_type(const Multi_parameter_generator<T>& f, std::size_t numberOfParameters) {
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
    std::vector<U> values(numberOfParameters, f[0]);
    for (std::size_t p = 0; p < f.size(); ++p) {
      values[p] = f[p];
    }
    return Out_multi_filtration(std::move(values), numberOfParameters);
  } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
    typename Out_multi_filtration::Underlying_container values(1, f.template as_type<U>());
    return Out_multi_filtration(std::move(values), numberOfParameters);
  } else if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    if (f.size() != 2 || numberOfParameters != 2)
      throw std::invalid_argument("Cannot convert a non-bifiltration to a bifiltration.");
    U inf = co ? Out_multi_filtration::T_m_inf : Out_multi_filtration::T_inf;
    std::vector<U> values(f[1] + 1, inf);
    values[f[1]] = f[0];
    return Out_multi_filtration(std::move(values), 2);
  } else {
    throw std::invalid_argument("Given out multi filtration value is not available.");
  }
}

/**
 * @ingroup multi_filtration
 *
 * @brief Converts the given multi filtration value into the type given as template argument.
 *
 * @tparam Out_multi_filtration New filtration value type. Has to be either
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration with desired template arguments.
 * @tparam T First template argument of the initial filtration value type.
 * @tparam Co Second template argument of the initial filtration value type.
 * @tparam Ensure1Criticality Third template argument of the initial filtration value type.
 * @param f Filtration value to convert.
 */
template <class Out_multi_filtration, typename T, bool Co, bool Ensure1Criticality>
inline Out_multi_filtration as_type(const Degree_rips_bifiltration<T, Co, Ensure1Criticality>& f)
{
  using U = typename Out_multi_filtration::value_type;
  constexpr bool co = Out_multi_filtration::has_negative_cones();
  constexpr bool one_crit = Out_multi_filtration::ensures_1_criticality();

  if constexpr (std::is_same_v<Out_multi_filtration, Degree_rips_bifiltration<U, co, one_crit> >) {
    return f.template as_type<U, co, one_crit>();
  } else {
    auto gen_index = [&f](std::size_t i) {
      if constexpr (Out_multi_filtration::has_negative_cones()) {
        return f.num_generators() - 1 - i;
      } else {
        return i;
      }
    };

    auto strictly_dominates = [](T a, T b) {
      if constexpr (Out_multi_filtration::has_negative_cones()) {
        return a < b;
      } else {
        return a > b;
      }
    };

    if (f.size() == 0) return Out_multi_filtration(0);

    std::vector<std::size_t> order(f.num_generators());
    std::iota(order.begin(), order.end(), 0);
    // lexicographical order
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j){
      if (f(i, 0) == f(j, 0)) return f(i, 1) < f(j, 1); // f(i, 1) and f(j, 1) cannot be equal for i != j
      return f(i, 0) < f(j, 0);
    });

    if constexpr (std::is_same_v<Out_multi_filtration, Multi_parameter_filtration<U, co, one_crit> >) {
      std::vector<U> values;
      values.reserve(f.num_generators() * 2);
      std::size_t g = order[gen_index(0)];
      T threshold = g;
      values.push_back(f(g, 0));
      values.push_back(threshold);
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        g = order[gen_index(i)];
        if (strictly_dominates(threshold, g)) {
          threshold = g;
          values.push_back(f(g, 0));
          values.push_back(threshold);
        }
      }
      if constexpr (co) {
        // lexicographical order
        const std::size_t max_idx = values.size() - 1;
        for (std::size_t i = 0; i < values.size() / 2; i += 2) {
          std::swap(values[i], values[max_idx - 1 - i]);
          std::swap(values[i + 1], values[max_idx - i]);
        }
      }

      return Out_multi_filtration(std::move(values), 2);
    } else if constexpr (std::is_same_v<Out_multi_filtration, Dynamic_multi_parameter_filtration<U, co, one_crit> >) {
      std::vector<Multi_parameter_generator<U> > values;
      values.reserve(f.num_generators());
      std::size_t g = order[gen_index(0)];
      T threshold = g;
      values.emplace_back(std::vector<T>{static_cast<T>(f(g, 0)), threshold});
      for (std::size_t i = 1; i < f.num_generators(); ++i) {
        g = order[gen_index(i)];
        if (strictly_dominates(threshold, g)) {
          threshold = g;
          std::vector<T> v = {static_cast<T>(f(g, 0)), threshold};
          values.emplace_back(std::move(v));
        }
      }
      if constexpr (co) {
        // lexicographical order
        std::reverse(values.begin(), values.end());
      }

      return Out_multi_filtration(std::move(values), 2);
    } else {
      throw std::invalid_argument("Given out multi filtration value is not available.");
    }
  }
}

/**
 * @private
 * Tries as good as it can to be precise, but there is probably room for improvements...
 */
template <typename T, typename U>
bool _safe_equal(T t, U u) {
  static_assert(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>, "_safe_equal requires arithmetic types");

  if constexpr (std::is_same_v<T, U>) {
    return t == u;
  } else if constexpr (std::is_integral_v<T> && std::is_integral_v<U>) {
    // integer x integer
    // in C++20, can be replaced with std::cmp_equal
    return boost::safe_numerics::safe<T>(t) == boost::safe_numerics::safe<U>(u);
  } else if constexpr (std::is_floating_point_v<T> && std::is_floating_point_v<U>) {
    // floating x floating
    using Wider = std::conditional_t<(sizeof(T) >= sizeof(U)), T, U>;
    return static_cast<Wider>(t) == static_cast<Wider>(u);
  } else {
    if constexpr (std::is_floating_point_v<T>) {
      // floating x integer
      // if the floating type is nan or has non-zero decimal values, it cannot be equal to an integer.
      if (std::isnan(t) || std::trunc(t) != t) return false;
      // the floating type does not have any decimal, i.e. is comparable with an integer
      if (t < static_cast<T>(std::numeric_limits<U>::min()) || t > static_cast<T>(std::numeric_limits<U>::max())) {
        return false;
      }
      return static_cast<U>(t) == u;
    } else {
      // integer x floating
      if (std::isnan(u) || std::trunc(u) != u) return false;
      if (u < static_cast<U>(std::numeric_limits<T>::min()) || u > static_cast<U>(std::numeric_limits<T>::max())) {
        return false;
      }
      return t == static_cast<T>(u);
    }
  }
}

/**
 * @brief Compares for equality two filtration values potentially not of the same type.
 * 
 * @tparam MultiFiltrationValue1 A filtration value type with methods: num_parameters(), num_generators()
 * and operator(g, p).
 * @tparam MultiFiltrationValue2 A filtration value type with methods: num_parameters(), num_generators()
 * and operator(g, p).
 */
template <class MultiFiltrationValue1, class MultiFiltrationValue2>
inline bool are_equal_filtration_values(const MultiFiltrationValue1& f1, const MultiFiltrationValue2& f2) {
  if constexpr (std::is_same_v<MultiFiltrationValue1, MultiFiltrationValue2>) {
    if (&f1 == &f2) return true;
  }
  if (f1.num_parameters() != f2.num_parameters()) return false;
  if (f1.num_generators() != f2.num_generators()) return false;

  for (std::size_t p = 0; p < f1.num_parameters(); ++p) {
    for (std::size_t g = 0; g < f1.num_generators(); ++g) {
      if (!_safe_equal(f1(g, p), f2(g, p))) return false;
    }
  }
  return true;
}

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_CONVERSIONS_H_
