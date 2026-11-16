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
 * @file Z2_field_element.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Z2_field_element class.
 */

#ifndef GUDHI_MATRIX_FIELD_Z2_H_
#define GUDHI_MATRIX_FIELD_Z2_H_

#include <utility>  // std::swap

#include <gudhi/Fields/Z2_field_operators.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Z2_field_element Z2_field_element.h gudhi/Fields/Z2_field_element.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of the \f$ \mathbb{F}_2 \f$ field.
 */
class Z2_field_element {
 public:
  using Op = Z2_field_operators;
  using Element = Op::Element;               /**< Type for the elements in the field. */
  using Characteristic = Op::Characteristic; /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor.
   */
  Z2_field_element() : element_(false) {}

  /**
   * @brief Constructor setting the element to the given value.
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Z2_field_element(Integer_type element) : element_(Op::get_value(element)) {}

  /**
   * @brief Copy constructor.
   *
   * @param toCopy Element to copy.
   */
  Z2_field_element(const Z2_field_element& toCopy) = default;

  /**
   * @brief Move constructor.
   *
   * @param toMove Element to move.
   */
  Z2_field_element(Z2_field_element&& toMove) noexcept : element_(std::exchange(toMove.element_, false)) {}

  ~Z2_field_element() = default;

  /**
   * @brief operator+=
   */
  friend void operator+=(Z2_field_element& f1, Z2_field_element const& f2) {
    Op::add_inplace(f1.element_, f2.element_);
  }

  /**
   * @brief operator+
   */
  friend Z2_field_element operator+(Z2_field_element f1, Z2_field_element const& f2) {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Z2_field_element& f, const Integer_type& v) {
    Op::add_inplace(f.element_, Op::get_value(v));
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator+(Z2_field_element f, const Integer_type& v) {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Z2_field_element f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Z2_field_element& f1, Z2_field_element const& f2) {
    Op::subtract_inplace_front(f1.element_, f2.element_);
  }

  /**
   * @brief operator-
   */
  friend Z2_field_element operator-(Z2_field_element f1, Z2_field_element const& f2) {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Z2_field_element& f, const Integer_type& v) {
    Op::subtract_inplace_front(f.element_, Op::get_value(v));
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator-(Z2_field_element f, const Integer_type& v) {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Z2_field_element& f) {
    return Op::subtract(Op::get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Z2_field_element& f1, Z2_field_element const& f2) {
    Op::multiply_inplace(f1.element_, f2.element_);
  }

  /**
   * @brief operator*
   */
  friend Z2_field_element operator*(Z2_field_element f1, Z2_field_element const& f2) {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Z2_field_element& f, const Integer_type& v) {
    Op::multiply_inplace(f.element_, Op::get_value(v));
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator*(Z2_field_element f, const Integer_type& v) {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Z2_field_element f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Z2_field_element& f1, const Z2_field_element& f2) { return f1.element_ == f2.element_; }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type v, const Z2_field_element& f) {
    return Op::get_value(v) == f.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Z2_field_element& f, const Integer_type v) {
    return Op::get_value(v) == f.element_;
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Z2_field_element& f1, const Z2_field_element& f2) { return !(f1 == f2); }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Z2_field_element& f) {
    return !(v == f);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Z2_field_element& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Z2_field_element& operator=(Z2_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  }

  /**
   * @brief Move assign operator.
   */
  Z2_field_element& operator=(Z2_field_element&& other) noexcept {
    if (this == &other) return *this;

    element_ = std::exchange(other.element_, false);
    return *this;
  }

  /**
   * @brief Assign operator.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Z2_field_element& operator=(Integer_type value) {
    element_ = Op::get_value(value);
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Z2_field_element& f1, Z2_field_element& f2) noexcept { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return static_cast<unsigned int>(element_); }

  /**
   * @brief Returns the opposite of the element.
   */
  [[nodiscard]] Z2_field_element get_opposite() const { return Op::get_opposite(element_); }

  /**
   * @brief Returns the inverse of the element.
   */
  [[nodiscard]] Z2_field_element get_inverse() const { return Op::get_inverse(element_); }

  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   *
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse and the second element is @p productOfCharacteristics.
   */
  [[nodiscard]] std::pair<Z2_field_element, Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    return Op::get_partial_inverse(element_, productOfCharacteristics);
  }

  /**
   * @brief Returns the additive identity of the field.
   *
   * @return false.
   */
  static Z2_field_element get_additive_identity() { return {Op::get_additive_identity()}; }

  /**
   * @brief Returns the multiplicative identity of the field.
   *
   * @return true.
   */
  static Z2_field_element get_multiplicative_identity() { return {Op::get_multiplicative_identity()}; }

  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   *
   * @param productOfCharacteristics Some value.
   * @return true.
   */
  static Z2_field_element get_partial_multiplicative_identity(Characteristic productOfCharacteristics) {
    return Op::get_partial_multiplicative_identity(productOfCharacteristics);
  }

  /**
   * @brief Returns the characteristic of the field, that is `2`.
   *
   * @return 2.
   */
  static constexpr Characteristic get_characteristic() { return Op::get_characteristic(); }

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  [[nodiscard]] Element get_value() const { return element_; }

 private:
  Element element_;
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_Z2_H_
