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
 * @file Multi_field_small_shared.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Shared_multi_field_element_with_small_characteristics class.
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_SHARED_H_
#define MATRIX_FIELD_MULTI_SMALL_SHARED_H_

#include <utility>  // std::swap

#include <gudhi/Debug_utils.h>
#include <gudhi/arithmetic.h>
#include <gudhi/Fields/Multi_field_small_operators.h>

namespace Gudhi {
namespace persistence_fields {

// TODO: either remove Unsigned_integer_type template or add it to Multi_field_operators_with_small_characteristics
/**
 * @class Shared_multi_field_element_with_small_characteristics Multi_field_small_shared.h \
 * gudhi/Fields/Multi_field_small_shared.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field, such that `productOfAllCharacteristics ^ 2` fits into
 * the given @p Unsigned_integer_type template argument. If each instantiation of the class can represent another
 * element, they all share the same characteristics. That is if the characteristics are set for one, they will be
 * set for all the others. The characteristics can be set before instantiating the elements with the static
 * @ref Shared_multi_field_element_with_small_characteristics::initialize method.
 *
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Shared_multi_field_element_with_small_characteristics {
 public:
  using Element = Multi_field_operators_with_small_characteristics::Element; /**< Type for the elements in the field. */
  using Characteristic =
      Multi_field_operators_with_small_characteristics::Characteristic; /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Shared_multi_field_element_with_small_characteristics() : element_(0) {}

  /**
   * @brief Constructor setting the element to the given value.
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if negative.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Shared_multi_field_element_with_small_characteristics(Integer_type element) : element_(op_.get_value(element)) {}

  /**
   * @brief Initialize the multi-field to the characteristics (primes) contained in the given interval.
   * Should be called first before constructing the field elements.
   * The characteristics must be small enough such that `productOfAllCharacteristics ^ 2` fits into an unsigned int.
   *
   * @param minimum Lowest value in the interval.
   * @param maximum Highest value in the interval.
   */
  static void initialize(int minimum, int maximum) { op_.set_characteristic(minimum, maximum); }

  /**
   * @brief operator+=
   */
  friend void operator+=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    op_.add_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator+
   */
  friend Shared_multi_field_element_with_small_characteristics operator+(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.add_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator+(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Shared_multi_field_element_with_small_characteristics f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    op_.subtract_inplace_front_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator-
   */
  friend Shared_multi_field_element_with_small_characteristics operator-(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.subtract_inplace_front_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator-(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Shared_multi_field_element_with_small_characteristics& f) {
    return op_.subtract_reduced(op_.get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    op_.multiply_inplace_reduced(f1.element_, f2.element_);
  }

  /**
   * @brief operator*
   */
  friend Shared_multi_field_element_with_small_characteristics operator*(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    op_.multiply_inplace_reduced(f.element_, op_.get_value(v));
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator*(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Shared_multi_field_element_with_small_characteristics f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f1,
                         const Shared_multi_field_element_with_small_characteristics& f2) {
    return f1.element_ == f2.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type v, const Shared_multi_field_element_with_small_characteristics& f) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator==
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return op_.get_value(v) == f.element_;
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f1,
                         const Shared_multi_field_element_with_small_characteristics& f2) {
    return !(f1 == f2);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Shared_multi_field_element_with_small_characteristics& f) {
    return !(v == f);
  }

  /**
   * @brief operator!=
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Shared_multi_field_element_with_small_characteristics& f1,
                   Shared_multi_field_element_with_small_characteristics& f2) noexcept {
    std::swap(f1.element_, f2.element_);
  }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the opposite of the element in the current product of all field.
   */
  Shared_multi_field_element_with_small_characteristics get_opposite() const { return op_.get_opposite(element_); }

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   */
  Shared_multi_field_element_with_small_characteristics get_inverse() const { return op_.get_inverse(element_); }

  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  std::pair<Shared_multi_field_element_with_small_characteristics, Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    return op_.get_partial_inverse(element_, productOfCharacteristics);
  }

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Shared_multi_field_element_with_small_characteristics get_additive_identity() {
    return {Multi_field_operators_with_small_characteristics::get_additive_identity()};
  }

  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Shared_multi_field_element_with_small_characteristics get_multiplicative_identity() {
    return {Multi_field_operators_with_small_characteristics::get_multiplicative_identity()};
  }

  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Shared_multi_field_element_with_small_characteristics get_partial_multiplicative_identity(
      const Characteristic& productOfCharacteristics) {
    return op_.get_partial_multiplicative_identity(productOfCharacteristics);
  }

  /**
   * @brief Returns the product of all characteristics.
   *
   * @return The product of all characteristics.
   */
  static Characteristic get_characteristic() { return op_.get_characteristic(); }

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  [[nodiscard]] Element get_value() const { return element_; }

 private:
  Element element_; /**< Element. */
  static inline Multi_field_operators_with_small_characteristics op_;
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_SMALL_SHARED_H_
