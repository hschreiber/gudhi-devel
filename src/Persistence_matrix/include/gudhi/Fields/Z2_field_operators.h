/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Z2_field_operators.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Z2_field_operators class.
 */

#ifndef GUDHI_MATRIX_FIELD_Z2_OPERATORS_H_
#define GUDHI_MATRIX_FIELD_Z2_OPERATORS_H_

#include <type_traits>
#include <utility>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Z2_field_operators Z2_field_operators.h gudhi/Fields/Z2_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the \f$ \mathbb{F}_2 \f$ field.
 */
class Z2_field_operators {
 public:
  using Element = bool;                /**< Type for the elements in the field. */
  using Characteristic = unsigned int; /**< Type for the field characteristic. */
  template <class T>
  using isUnsignedInteger = std::enable_if_t<std::is_unsigned_v<T> >;
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;
  template <typename... Ts>
  using isIntegers = std::enable_if_t<(std::is_integral_v<Ts> && ...)>;

  constexpr static const Characteristic nullCharacteristic = 2;

  /** @brief Returns the characteristic of the field, that is `2`. */
  static constexpr Characteristic get_characteristic() noexcept { return 2; }

  /**
   * @brief Returns the value of an integer in the field, that is its parity.
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   * @param e Integer to return the value from.
   * @return A boolean representing `e % 2` (non negative, also for negative @p e).
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  static constexpr Element get_value(Integer_type e) noexcept {
    if constexpr (std::is_same_v<Integer_type, bool>) {
      return e;
    } else if constexpr (std::is_signed_v<Integer_type>) {
      return e % 2 != 0;
    } else {
      return (e & 1U) != 0;
    }
  }

  /** @brief Returns `(e1 + e2) % 2`. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr Element add(I1 e1, I2 e2) noexcept {
    return get_value(e1) != get_value(e2);
  }

  /** @brief Stores `(e1 + e2) % 2` in @p e1. */
  template <typename I, class = isIntegers<I> >
  static constexpr void add_inplace(Element& e1, I e2) noexcept {
    e1 = e1 != get_value(e2);
  }

  /** @brief Returns `(e1 - e2) % 2`. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr Element subtract(I1 e1, I2 e2) noexcept {
    return get_value(e1) != get_value(e2);
  }

  /** @brief Stores `(e1 - e2) % 2` in @p e1. */
  template <typename I, class = isIntegers<I> >
  static constexpr void subtract_inplace_front(Element& e1, I e2) noexcept {
    e1 = e1 != get_value(e2);
  }

  /** @brief Stores `(e1 - e2) % 2` in @p e2. */
  template <typename I, class = isIntegers<I> >
  static constexpr void subtract_inplace_back(I e1, Element& e2) noexcept {
    e2 = get_value(e1) != e2;
  }

  /** @brief Returns `(e1 * e2) % 2`. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr Element multiply(I1 e1, I2 e2) noexcept {
    return get_value(e1) && get_value(e2);
  }

  /** @brief Stores `(e1 * e2) % 2` in @p e1. */
  template <typename I, class = isIntegers<I> >
  static constexpr void multiply_inplace(Element& e1, I e2) noexcept {
    e1 = e1 && get_value(e2);
  }

  /** @brief Returns `(e * m + a) % 2`. */
  template <typename I1, typename I2, typename I3, class = isIntegers<I1, I2, I3> >
  static constexpr Element multiply_and_add(I1 e, I2 m, I3 a) noexcept {
    return (get_value(e) && get_value(m)) != get_value(a);
  }

  /** @brief Stores `(e * m + a) % 2` in @p e. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr void multiply_and_add_inplace_front(Element& e, I1 m, I2 a) noexcept {
    e = (e && get_value(m)) != get_value(a);
  }

  /** @brief Stores `(e * m + a) % 2` in @p a. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr void multiply_and_add_inplace_back(I1 e, I2 m, Element& a) noexcept {
    a = (get_value(e) && get_value(m)) != a;
  }

  /** @brief Returns `((e + a) * m) % 2`. */
  template <typename I1, typename I2, typename I3, class = isIntegers<I1, I2, I3> >
  static constexpr Element add_and_multiply(I1 e, I2 a, I3 m) noexcept {
    return (get_value(e) != get_value(a)) && get_value(m);
  }

  /** @brief Stores `((e + a) * m) % 2` in @p e. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr void add_and_multiply_inplace_front(Element& e, I1 a, I2 m) noexcept {
    e = (e != get_value(a)) && get_value(m);
  }

  /** @brief Stores `((e + a) * m) % 2` in @p m. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr void add_and_multiply_inplace_back(I1 e, I2 a, Element& m) noexcept {
    m = (get_value(e) != get_value(a)) && m;
  }

  /** @brief Returns true if `e1 % 2 == e2 % 2`. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  static constexpr bool are_equal(I1 e1, I2 e2) noexcept {
    return get_value(e1) == get_value(e2);
  }

  /** @brief Returns the opposite of `e % 2`, that is `e % 2`. */
  template <typename I, class = isIntegers<I> >
  static constexpr Element get_opposite(I e) noexcept {
    return get_value(e);
  }

  /** @brief Replaces @p e by its opposite, that is itself. */
  static constexpr void get_opposite_inplace([[maybe_unused]] Element& e) noexcept {}

  /** @brief Returns the inverse of `e % 2`, that is `e % 2`. (The inverse of 0 is 0 by convention.) */
  template <typename I, class = isIntegers<I> >
  static constexpr Element get_inverse(I e) noexcept {
    return get_value(e);
  }

  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   */
  template <typename I, class = isIntegers<I> >
  static constexpr std::pair<Element, Characteristic> get_partial_inverse(
      I e, Characteristic productOfCharacteristics) noexcept {
    return {get_inverse(e), productOfCharacteristics};
  }

  /** @brief Returns the additive identity of the field, that is 0. */
  static constexpr Element get_additive_identity() noexcept { return false; }

  /** @brief Returns the multiplicative identity of the field, that is 1. */
  static constexpr Element get_multiplicative_identity() noexcept { return true; }

  /** @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field. */
  static constexpr Element get_partial_multiplicative_identity([[maybe_unused]] Characteristic p) noexcept {
    return true;
  }

  // ---------------------------------------------------------------------------------------------------------
  // "_reduced" variants: for interface uniformity with the other fields (where they skip the reduction of the
  // arguments). Here, elements are booleans and are always reduced, so they are plain aliases.

  /** @brief Same as @ref add, but arguments are already of type @ref Element and reduced. */
  static constexpr Element add_reduced(Element e1, Element e2) noexcept { return e1 != e2; }

  /** @brief Same as @ref add_inplace, but arguments are already of type @ref Element and reduced. */
  static constexpr void add_inplace_reduced(Element& e1, Element e2) noexcept { e1 = e1 != e2; }

  /** @brief Same as @ref subtract, but arguments are already of type @ref Element and reduced. */
  static constexpr Element subtract_reduced(Element e1, Element e2) noexcept { return e1 != e2; }

  /** @brief Same as @ref subtract_inplace, but arguments are already of type @ref Element and reduced. */
  static constexpr void subtract_inplace_front_reduced(Element& e1, Element e2) noexcept { e1 = e1 != e2; }

  /** @brief Same as @ref multiply, but arguments are already of type @ref Element and reduced. */
  static constexpr Element multiply_reduced(Element e1, Element e2) noexcept { return e1 && e2; }

  /** @brief Same as @ref multiply_inplace, but arguments are already of type @ref Element and reduced. */
  static constexpr void multiply_inplace_reduced(Element& e1, Element e2) noexcept { e1 = e1 && e2; }

  /** @brief Same as @ref multiply_and_add, but arguments are already of type @ref Element and reduced. */
  static constexpr Element multiply_and_add_reduced(Element e, Element m, Element a) noexcept { return (e && m) != a; }

  /** @brief Same as @ref multiply_and_add_inplace_front, but arguments are already of type @ref Element and reduced. */
  static constexpr void multiply_and_add_inplace_front_reduced(Element& e, Element m, Element a) noexcept {
    e = (e && m) != a;
  }

  /** @brief Same as @ref multiply_and_add_inplace_front, but arguments are already of type @ref Element and reduced. */
  static constexpr void multiply_and_add_inplace_back_reduced(Element e, Element m, Element& a) noexcept {
    a = (e && m) != a;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_Z2_OPERATORS_H_
