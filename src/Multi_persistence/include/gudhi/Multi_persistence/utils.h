/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @private
 * @file utils.h
 * @author Hannah Schreiber
 */

#ifndef MP_UTILS_H_
#define MP_UTILS_H_

#include <iterator>
#include <type_traits>

namespace Gudhi {

namespace multi_persistence {

namespace details {

template <class T, typename = void>
struct is_forward_iterator : std::false_type {};

template <class T>
struct is_forward_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category>>
    : std::bool_constant<
          std::is_base_of_v<std::forward_iterator_tag, typename std::iterator_traits<T>::iterator_category>> {};

}  // namespace details

template <class T>
constexpr bool is_forward_iterator_v = details::is_forward_iterator<T>::value;

template <typename T>
struct type_identity {
  using type = T;
};

// std::make_signed_t does not compile for T signed and std::conditional evaluates both possibilities
// so this trick is necessary if we want to avoid using `if constexpr` everywhere
template <typename T>
using maybe_make_signed_t =
    typename std::conditional_t<std::is_unsigned_v<T>, std::make_signed<T>, type_identity<T>>::type;

}  // namespace multi_persistence

}  // namespace Gudhi

#endif  // MP_UTILS_H_
