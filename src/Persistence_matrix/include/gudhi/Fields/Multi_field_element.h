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
 * @file Multi_field_element.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_element class.
 */

#ifndef GUDHI_MATRIX_FIELD_MULTI_H_
#define GUDHI_MATRIX_FIELD_MULTI_H_

#include <utility>  // std::swap

#include <gudhi/Debug_utils.h>
#include <gudhi/arithmetic.h>
#include <gudhi/Fields/Multi_field_operators.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_element Multi_field_element.h gudhi/Fields/Multi_field_element.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field.
 * The characteristics will corresponds to all prime numbers in the interval given as template.
 *
 * @tparam minimum Interval closed lower bound.
 * @tparam maximum Interval closed upper bound.
 */
template <unsigned int minimum, unsigned int maximum>
class Multi_field_element {
 public:
  using Element = Multi_field_operators::Element;               /**< Type for the elements in the field. */
  using Characteristic = Multi_field_operators::Characteristic; /**< Type for the field characteristic. */

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Multi_field_element() : element_(0) {
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
  Multi_field_element(const Element& element) : element_(op_.get_value(element)) {
    static_assert(maximum >= 2, "Characteristics has to be positive.");
    static_assert(minimum <= maximum, "The given interval is not valid.");
    static_assert(minimum != maximum || Gudhi::is_prime(minimum),
                  "The given interval does not contain a prime number.");
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element& f1, Multi_field_element const& f2) {
    op_.add_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator+
   */
  friend Multi_field_element operator+(Multi_field_element f1, Multi_field_element const& f2) {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element& f, const Element& v) {
    op_.add_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator+
   */
  friend Multi_field_element operator+(Multi_field_element f, const Element& v) {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   */
  friend Element operator+(const Element& v, Multi_field_element f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element& f1, Multi_field_element const& f2) {
    op_.subtract_inplace_front_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator-
   */
  friend Multi_field_element operator-(Multi_field_element f1, Multi_field_element const& f2) {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element& f, const Element& v) {
    op_.subtract_inplace_front_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator-
   */
  friend Multi_field_element operator-(Multi_field_element f, const Element& v) {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   */
  friend Element operator-(const Element& v, Multi_field_element const& f) {
    return op_.subtract_reduced(op_.get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element& f1, Multi_field_element const& f2) {
    op_.multiply_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator*
   */
  friend Multi_field_element operator*(Multi_field_element f1, Multi_field_element const& f2) {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element& f, const Element& v) {
    op_.multiply_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator*
   */
  friend Multi_field_element operator*(Multi_field_element f, const Element& v) {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   */
  friend Element operator*(const Element& v, Multi_field_element f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element& f1, const Multi_field_element& f2) {
    return f1.element_ == f2.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Element& v, const Multi_field_element& f) { return op_.get_value(v) == f.element_; }

  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element& f, const Element& v) { return op_.get_value(v) == f.element_; }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element& f1, const Multi_field_element& f2) { return !(f1 == f2); }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Element& v, const Multi_field_element& f) { return !(v == f); }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element& f, const Element& v) { return !(v == f); }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_element& f1, Multi_field_element& f2) noexcept { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_.get_ui(); }

  /**
   * @brief Casts the element into an mpz_class.
   */
  explicit operator mpz_class() const { return element_; }

  /**
   * @brief Returns the opposite of the element in the current product of all field.
   */
  Multi_field_element get_opposite() const { return op_.get_opposite(element_); }

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   */
  [[nodiscard]] Multi_field_element get_inverse() const { return op_.get_inverse(element_); }

  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  [[nodiscard]] std::pair<Multi_field_element, Characteristic> get_partial_inverse(
      const Characteristic& productOfCharacteristics) const {
    return op_.get_partial_inverse(element_, productOfCharacteristics);
  }

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Multi_field_element get_additive_identity() { return {Multi_field_operators::get_additive_identity()}; }

  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Multi_field_element get_multiplicative_identity() {
    return {Multi_field_operators::get_multiplicative_identity()};
  }

  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Multi_field_element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics) {
    return op_.get_partial_multiplicative_identity(productOfCharacteristics);
  }

  /**
   * @brief Returns the product of all characteristics.
   */
  static Characteristic get_characteristic() { return op_.get_characteristic(); }

  /**
   * @brief Returns the value of the element.
   */
  [[nodiscard]] Element get_value() const { return element_; }

 private:
  Element element_;
  static inline const Multi_field_operators op_ = {minimum, maximum};
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_MULTI_H_
