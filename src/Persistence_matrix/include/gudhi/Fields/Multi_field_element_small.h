/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_field_element_small.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_element_with_small_characteristics class.
 */

#ifndef GUDHI_MATRIX_FIELD_MULTI_SMALL_H_
#define GUDHI_MATRIX_FIELD_MULTI_SMALL_H_

#include <utility>  // std::swap

#include <gudhi/Debug_utils.h>
#include <gudhi/arithmetic.h>
#include <gudhi/Fields/Multi_field_small_operators.h>

namespace Gudhi {
namespace persistence_fields {

// TODO: either remove Unsigned_integer_type template or add it to Multi_field_operators_with_small_characteristics
/**
 * @class Multi_field_element_with_small_characteristics Multi_field_element_small.h
 * gudhi/Fields/Multi_field_element_small.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field, such that the product of all characteristics fits into
 * the given @p Unsigned_integer_type template argument. The characteristics will corresponds to all prime numbers
 * in the interval given as other template arguments.
 *
 * @tparam minimum Interval closed lower bound.
 * @tparam maximum Interval closed upper bound.
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <unsigned int minimum, unsigned int maximum, typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Multi_field_element_with_small_characteristics {
 public:
  using Element = Multi_field_operators_with_small_characteristics::Element; /**< Type for the elements in the field. */
  using Characteristic =
      Multi_field_operators_with_small_characteristics::Characteristic; /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Multi_field_element_with_small_characteristics() : element_(0) {
    static_assert(maximum >= 2, "Characteristics have to be positive.");
    static_assert(minimum <= maximum, "The given interval is not valid.");
    static_assert(minimum != maximum || Gudhi::is_prime(minimum),
                  "The given interval does not contain a prime number.");
  }

  /**
   * @brief Constructor setting the element to the given value.
   *
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Multi_field_element_with_small_characteristics(Integer_type element) : element_(op_.get_value(element)) {
    static_assert(maximum >= 2, "Characteristics has to be positive.");
    static_assert(minimum <= maximum, "The given interval is not valid.");
    static_assert(minimum != maximum || Gudhi::is_prime(minimum),
                  "The given interval does not contain a prime number.");
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    op_.add_inplace(f1.element_, f2.element_);
  }

  /**
   * @brief operator+
   */
  friend Multi_field_element_with_small_characteristics operator+(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.add_inplace(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator+(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Multi_field_element_with_small_characteristics f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    op_.subtract_inplace_front(f1.element_, f2.element_);
  }

  /**
   * @brief operator-
   */
  friend Multi_field_element_with_small_characteristics operator-(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.subtract_inplace_front(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator-(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Multi_field_element_with_small_characteristics& f) {
    return op_.subtract(op_.get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    op_.multiply_inplace(f1.element_, f2.element_);
  }

  /**
   * @brief operator*
   */
  friend Multi_field_element_with_small_characteristics operator*(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.multiply_inplace(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator*(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Multi_field_element_with_small_characteristics f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element_with_small_characteristics& f1,
                         const Multi_field_element_with_small_characteristics& f2) {
    return f1.element_ == f2.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type v, const Multi_field_element_with_small_characteristics& f) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element_with_small_characteristics& f1,
                         const Multi_field_element_with_small_characteristics& f2) {
    return !(f1 == f2);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Multi_field_element_with_small_characteristics& f) {
    return !(v == f);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_element_with_small_characteristics& f1,
                   Multi_field_element_with_small_characteristics& f2) noexcept {
    std::swap(f1.element_, f2.element_);
  }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the opposite of the element in the current product of all field.
   */
  Multi_field_element_with_small_characteristics get_opposite() const { return op_.get_opposite(element_); }

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   */
  Multi_field_element_with_small_characteristics get_inverse() const { return op_.get_inverse(element_); }

  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  std::pair<Multi_field_element_with_small_characteristics, Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    return op_.get_partial_inverse(element_, productOfCharacteristics);
  }

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Multi_field_element_with_small_characteristics get_additive_identity() {
    return {Multi_field_operators_with_small_characteristics::get_additive_identity()};
  }

  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Multi_field_element_with_small_characteristics get_multiplicative_identity() {
    return {Multi_field_operators_with_small_characteristics::get_multiplicative_identity()};
  }

  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Multi_field_element_with_small_characteristics get_partial_multiplicative_identity(
      Characteristic productOfCharacteristics) {
    return op_.get_partial_multiplicative_identity(productOfCharacteristics);
  }

  /**
   * @brief Returns the product of all characteristics.
   *
   * @return The product of all characteristics.
   */
  static constexpr Characteristic get_characteristic() { return op_.get_characteristic(); }

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  [[nodiscard]] Element get_value() const { return element_; }

 private:
  Element element_;
  static inline const Multi_field_operators_with_small_characteristics op_ = {minimum, maximum};
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_MULTI_SMALL_H_
