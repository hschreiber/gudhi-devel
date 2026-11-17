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
 * @file Zp_field_element.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Zp_field_element class.
 */

#ifndef GUDHI_MATRIX_FIELD_ZP_H_
#define GUDHI_MATRIX_FIELD_ZP_H_

#include <utility>

#include <gudhi/arithmetic.h>
#include <gudhi/Fields/Zp_field_operators.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Zp_field_element Zp_field_element.h gudhi/Fields/Zp_field_element.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 *
 * @tparam characteristic Value of the characteristic of the field. Has to be a positive prime number.
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <typename Unsigned_integer_type, Unsigned_integer_type characteristic,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Zp_field_element {
 public:
  using Op = Zp_field_operators<Unsigned_integer_type>;
  using Element = typename Op::Element;               /**< Type for the elements in the field. */
  using Characteristic = typename Op::Characteristic; /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Zp_field_element() : element_(0) {
    static_assert(Gudhi::is_prime(characteristic), "Characteristic has to be a prime number.");
  }

  /**
   * @brief Constructor setting the element to the given value.
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Zp_field_element(Integer_type element) : element_(op_.get_value(element)) {
    static_assert(Gudhi::is_prime(characteristic), "Characteristic has to be a prime number.");
  }

  /**
   * @brief Copy constructor.
   *
   * @param toCopy Element to copy.
   */
  Zp_field_element(const Zp_field_element& toCopy) = default;

  /**
   * @brief Move constructor.
   *
   * @param toMove Element to move.
   */
  Zp_field_element(Zp_field_element&& toMove) noexcept : element_(std::exchange(toMove.element_, 0)) {}

  ~Zp_field_element() = default;

  /**
   * @brief operator+=
   */
  friend void operator+=(Zp_field_element& f1, Zp_field_element const& f2) {
    op_.add_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator+
   */
  friend Zp_field_element operator+(Zp_field_element f1, Zp_field_element const& f2) {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Zp_field_element& f, const Integer_type& v) {
    op_.add_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator+(Zp_field_element f, const Integer_type& v) {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Zp_field_element f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Zp_field_element& f1, Zp_field_element const& f2) {
    op_.subtract_inplace_front_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator-
   */
  friend Zp_field_element operator-(Zp_field_element f1, Zp_field_element const& f2) {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Zp_field_element& f, const Integer_type& v) {
    op_.subtract_inplace_front_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator-(Zp_field_element f, const Integer_type& v) {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Zp_field_element& f) {
    return op_.subtract_reduced(op_.get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Zp_field_element& f1, Zp_field_element const& f2) {
    op_.multiply_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator*
   */
  friend Zp_field_element operator*(Zp_field_element f1, Zp_field_element const& f2) {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Zp_field_element& f, const Integer_type& v) {
    op_.multiply_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator*(Zp_field_element f, const Integer_type& v) {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Zp_field_element f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Zp_field_element& f1, const Zp_field_element& f2) { return f1.element_ == f2.element_; }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type v, const Zp_field_element& f) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Zp_field_element& f, const Integer_type v) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Zp_field_element& f1, const Zp_field_element& f2) { return !(f1 == f2); }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Zp_field_element& f) {
    return !(v == f);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Zp_field_element& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Zp_field_element& operator=(Zp_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  }

  /**
   * @brief Move assign operator.
   */
  Zp_field_element& operator=(Zp_field_element&& other) noexcept {
    if (this == &other) return *this;

    element_ = std::exchange(other.element_, 0);
    return *this;
  }

  /**
   * @brief Assign operator.
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Zp_field_element& operator=(const Integer_type& value) {
    element_ = op_.get_value(value);
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Zp_field_element& f1, Zp_field_element& f2) noexcept { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the opposite of the element.
   */
  [[nodiscard]] Zp_field_element get_opposite() const { return op_.get_opposite_reduced(element_); }

  /**
   * @brief Returns the inverse of the element.
   */
  [[nodiscard]] Zp_field_element get_inverse() const { return op_.get_inverse_reduced(element_); }

  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   *
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse and the second element is @p productOfCharacteristics.
   */
  [[nodiscard]] std::pair<Zp_field_element, Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    return op_.get_partial_inverse_reduced(element_, productOfCharacteristics);
  }

  /**
   * @brief Returns the additive identity of the field.
   *
   * @return 0.
   */
  static Zp_field_element get_additive_identity() { return {Op::get_additive_identity()}; }

  /**
   * @brief Returns the multiplicative identity of the field.
   *
   * @return 1.
   */
  static Zp_field_element get_multiplicative_identity() { return {Op::get_multiplicative_identity()}; }

  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   *
   * @param productOfCharacteristics Some value.
   * @return 1.
   */
  static Zp_field_element get_partial_multiplicative_identity(Characteristic productOfCharacteristics) {
    return Op::get_partial_multiplicative_identity(productOfCharacteristics);
  }

  /**
   * @brief Returns the current characteristic.
   *
   * @return The value of the current characteristic.
   */
  static Characteristic get_characteristic() { return op_.get_characteristic(); }

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  Element get_value() const { return element_; }

 private:
  Element element_; /**< Field element. */
  static inline const Op op_ = {characteristic};
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_ZP_H_
