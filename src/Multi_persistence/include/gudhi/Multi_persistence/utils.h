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

#include <cstddef>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <vector>

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

template <typename T, class>
inline char *serialize_value_to_char_buffer(T value, char *start);
template <typename T>
inline char *serialize_value_to_char_buffer(const std::vector<T> &range, char *start);
template <typename T1, typename T2>
inline char *serialize_value_to_char_buffer(const std::pair<T1, T2> &pair, char *start);

template <typename T, class>
inline const char *deserialize_value_from_char_buffer(T &value, const char *start);
template <typename T>
inline const char *deserialize_value_from_char_buffer(std::vector<T> &range, const char *start);
template <typename T1, typename T2>
inline const char *deserialize_value_from_char_buffer(std::pair<T1, T2> &pair, const char *start);

template <typename T, class>
inline std::size_t get_serialization_size_of([[maybe_unused]] T value);
template <typename T>
inline std::size_t get_serialization_size_of(const std::vector<T> &range);
template <typename T1, typename T2>
inline std::size_t get_serialization_size_of(const std::pair<T1, T2> &pair);

template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline char *serialize_value_to_char_buffer(T value, char *start) {
  const std::size_t argSize = sizeof(T);
  memcpy(start, &value, argSize);
  return start + argSize;
}

template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline const char *deserialize_value_from_char_buffer(T &value, const char *start) {
  const std::size_t argSize = sizeof(T);
  memcpy(&value, start, argSize);
  return start + argSize;
}

template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline std::size_t get_serialization_size_of([[maybe_unused]] T value) {
  return sizeof(T);
}

template <typename T>
inline char *serialize_value_to_char_buffer(const std::vector<T> &range, char *start) {
  const std::size_t sizeSize = sizeof(std::size_t);
  const std::size_t length = range.size();
  char *curr = start;
  memcpy(curr, &length, sizeSize);
  curr += sizeSize;
  for (const T &v : range) {
    curr = serialize_value_to_char_buffer(v, curr);
  }
  return curr;
}

template <typename T>
inline const char *deserialize_value_from_char_buffer(std::vector<T> &range, const char *start) {
  const std::size_t sizeSize = sizeof(std::size_t);
  std::size_t length;
  const char *curr = start;
  memcpy(&length, curr, sizeSize);
  curr += sizeSize;
  range.resize(length);
  for (T &v : range) {
    curr = deserialize_value_from_char_buffer(v, curr);
  }
  return curr;
}

template <typename T>
inline std::size_t get_serialization_size_of(const std::vector<T> &range) {
  std::size_t size = sizeof(std::size_t);
  for (const T &v : range) {
    size += get_serialization_size_of(v);
  }
  return size;
}

template <typename T1, typename T2>
inline char *serialize_value_to_char_buffer(const std::pair<T1, T2> &pair, char *start) {
  char *curr = start;
  curr = serialize_value_to_char_buffer(pair.first, curr);
  curr = serialize_value_to_char_buffer(pair.second, curr);
  return curr;
}

template <typename T1, typename T2>
inline const char *deserialize_value_from_char_buffer(std::pair<T1, T2> &pair, const char *start) {
  const char *curr = start;
  curr = deserialize_value_from_char_buffer(pair.first, curr);
  curr = deserialize_value_from_char_buffer(pair.second, curr);
  return curr;
}

template <typename T1, typename T2>
inline std::size_t get_serialization_size_of(const std::pair<T1, T2> &pair) {
  std::size_t size = get_serialization_size_of(pair.first);
  size += get_serialization_size_of(pair.second);
  return size;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_UTILS_H_
