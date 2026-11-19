/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Zp_field_operators.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Zp_field_operators class.
 */

#ifndef GUDHI_MATRIX_FIELD_ZP_OPERATOR_H_
#define GUDHI_MATRIX_FIELD_ZP_OPERATOR_H_

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <gudhi/arithmetic.h>

// Hint for the (very common) case of already reduced arguments, so that the compiler keeps a predictable branch
// instead of speculatively computing the reduction.
#if defined(__GNUC__) || defined(__clang__)
#define GUDHI_ZP_LIKELY(x) __builtin_expect(!!(x), 1)
#else
#define GUDHI_ZP_LIKELY(x) (x)
#endif

namespace Gudhi {
/// Field namespace
namespace persistence_fields {

/**
 * @class Zp_field_operators Zp_field_operators.h gudhi/Fields/Zp_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 *
 * Design notes:
 * - Elements are always stored fully reduced, that is in \f$ [0, p-1] \f$. Methods without suffix accept any native
 *   integer type (signed or not, up to 64 bits, a different type for each argument) and any value, also negative
 *   ones, and reduce their arguments first. They return an @ref Element. The "in place" versions store the result in
 *   an @ref Element reference (the value it contains before the call can be any @ref Element, reduced or not).
 *   Methods with the suffix `_reduced` only accept @ref Element arguments which are already in
 *   \f$ [0, p-1] \f$ (no check) and skip the reduction: use them in hot loops where this is an invariant.
 * - Nothing overflows: products are computed on the full width (see @ref Gudhi::numbers::mul_wide) and reduced with
 *   @ref Gudhi::numbers::Modular_reducer, which uses neither hardware division nor compiler support for 128 bits
 *   integers: Barrett reduction for 32 bits elements and division by an invariant integer for 64 bits elements.
 * - For repeated multiplications by the same element, see @ref prepare_multiplier.
 *
 * @tparam Unsigned_integer_type Either a 32 bits or a 64 bits native unsigned integer type (`bool` excluded).
 * Will be used as the field element type. Any other type is rejected at compile time. The characteristic has to be
 * strictly smaller than \f$ 2^{31} \f$ (resp. \f$ 2^{63} \f$) for 32 bits (resp. 64 bits) types.
 */
template <typename Unsigned_integer_type = unsigned int>
class Zp_field_operators {
  static_assert(std::is_unsigned_v<Unsigned_integer_type> && !std::is_same_v<Unsigned_integer_type, bool> &&
                    (sizeof(Unsigned_integer_type) == 4 || sizeof(Unsigned_integer_type) == 8),
                "Zp_field_operators supports only 32 bits and 64 bits unsigned integer types.");

 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;        /**< Type for the field characteristic. */
  template <typename... Ts>
  using isIntegers = std::enable_if_t<(std::is_integral_v<Ts> && ...)>;

  /**
   * @brief Value of a non initialized characteristic.
   */
  constexpr static const Characteristic nullCharacteristic = 0;

  /**
   * @brief Default constructor. If a non-zero characteristic is given, initializes the field with it.
   * The characteristic can later be changed again or initialized with @ref set_characteristic.
   *
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  Zp_field_operators(Characteristic characteristic = nullCharacteristic) : characteristic_(nullCharacteristic) {
    if (characteristic != nullCharacteristic) set_characteristic(characteristic);
  }

  /**
   * @brief Sets the characteristic of the field.
   *
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   * @throw std::invalid_argument If @p characteristic is not a prime number or is too large for @ref Element.
   * If the exception is thrown, the field is left unchanged.
   */
  void set_characteristic(Characteristic characteristic) {
    if (characteristic <= 1) throw std::invalid_argument("Characteristic must be a strictly positive prime number.");
    if (characteristic >= (Element(1) << (kBits - 1)))
      throw std::invalid_argument("Characteristic must be smaller than 2^(number of bits of the element type - 1).");
    if (!Gudhi::is_prime(characteristic)) throw std::invalid_argument("Characteristic must be a prime number.");

    std::vector<Element> inverse;
    if (characteristic <= kMaxTableCharacteristic) {
      // inv[i] = -(p / i) * inv[p % i]  (mod p), valid because p is prime.
      inverse.resize(characteristic);
      inverse[0] = 0;
      inverse[1] = 1;
      for (Element i = 2; i < characteristic; ++i) {
        const auto v = static_cast<Element>(
            (static_cast<std::uint64_t>(characteristic / i) * inverse[characteristic % i]) % characteristic);
        inverse[i] = characteristic - v;
      }
    }

    // nothing can throw anymore
    characteristic_ = characteristic;
    reducer_ = Reducer(characteristic);
    inverse_.swap(inverse);
  }

  /**
   * @brief Returns the current characteristic.
   */
  [[nodiscard]] const Characteristic& get_characteristic() const { return characteristic_; }

  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   *
   * @tparam Integer_type Any native integer type, signed or not, of any size.
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  template <typename Integer_type, class = std::enable_if_t<std::is_integral_v<Integer_type> > >
  [[nodiscard]] Element get_value(Integer_type e) const {
    if constexpr (std::is_signed_v<Integer_type>) {
      using U = std::make_unsigned_t<Integer_type>;
      if (e >= 0) return _reduce_unsigned(static_cast<U>(e));
      const Element r = _reduce_unsigned(static_cast<U>(U(0) - static_cast<U>(e)));  // reduces |e|
      return r == 0 ? Element(0) : characteristic_ - r;
    } else {
      return _reduce_unsigned(e);
    }
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations accepting arbitrary arguments.

  /** @brief Returns `(e1 + e2) % characteristic`, such that the result is positive. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  [[nodiscard]] Element add(I1 e1, I2 e2) const {
    return add_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Stores in the first element `(e1 + e2) % characteristic`, such that the result is positive. */
  template <typename I, class = isIntegers<I> >
  void add_inplace(Element& e1, I e2) const {
    e1 = add_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Returns `(e1 - e2) % characteristic`, such that the result is positive. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  [[nodiscard]] Element subtract(I1 e1, I2 e2) const {
    return subtract_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Stores in the first element `(e1 - e2) % characteristic`, such that the result is positive. */
  template <typename I, class = isIntegers<I> >
  void subtract_inplace_front(Element& e1, I e2) const {
    e1 = subtract_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Stores in the second element `(e1 - e2) % characteristic`, such that the result is positive. */
  template <typename I, class = isIntegers<I> >
  void subtract_inplace_back(I e1, Element& e2) const {
    e2 = subtract_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Returns `(e1 * e2) % characteristic`, such that the result is positive. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  [[nodiscard]] Element multiply(I1 e1, I2 e2) const {
    return multiply_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Stores in the first element `(e1 * e2) % characteristic`, such that the result is positive. */
  template <typename I, class = isIntegers<I> >
  void multiply_inplace(Element& e1, I e2) const {
    e1 = multiply_reduced(get_value(e1), get_value(e2));
  }

  /** @brief Returns `(e * m + a) % characteristic`, such that the result is positive. */
  template <typename I1, typename I2, typename I3, class = isIntegers<I1, I2, I3> >
  [[nodiscard]] Element multiply_and_add(I1 e, I2 m, I3 a) const {
    return multiply_and_add_reduced(get_value(e), get_value(m), get_value(a));
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p e. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  void multiply_and_add_inplace_front(Element& e, I1 m, I2 a) const {
    e = multiply_and_add_reduced(get_value(e), get_value(m), get_value(a));
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p a. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  void multiply_and_add_inplace_back(I1 e, I2 m, Element& a) const {
    a = multiply_and_add_reduced(get_value(e), get_value(m), get_value(a));
  }

  /** @brief Returns `((e + a) * m) % characteristic`, such that the result is positive. */
  template <typename I1, typename I2, typename I3, class = isIntegers<I1, I2, I3> >
  [[nodiscard]] Element add_and_multiply(I1 e, I2 a, I3 m) const {
    return multiply_reduced(add_reduced(get_value(e), get_value(a)), get_value(m));
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p e. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  void add_and_multiply_inplace_front(Element& e, I1 a, I2 m) const {
    e = multiply_reduced(add_reduced(get_value(e), get_value(a)), get_value(m));
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p m. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  void add_and_multiply_inplace_back(I1 e, I2 a, Element& m) const {
    m = multiply_reduced(add_reduced(get_value(e), get_value(a)), get_value(m));
  }

  /** @brief Returns true if `e1 % characteristic == e2 % characteristic`. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  [[nodiscard]] bool are_equal(I1 e1, I2 e2) const {
    return get_value(e1) == get_value(e2);
  }

  /** @brief Returns the opposite of `e % characteristic` in the field. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Element get_opposite(I e) const {
    const Element v = get_value(e);
    return v == 0 ? Element(0) : characteristic_ - v;
  }

  /** @brief Replaces @p e by its opposite in the field. */
  void get_opposite_inplace(Element& e) const { e = get_opposite(e); }

  /** @brief Returns the inverse of `e % characteristic` in the field. The inverse of 0 is 0 by convention. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Element get_inverse(I e) const {
    return get_inverse_reduced(get_value(e));
  }

  /** @brief For interface purposes with multi-fields. Returns the inverse together with the second argument. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse(I e,
                                                                       Characteristic productOfCharacteristics) const {
    return {get_inverse(e), productOfCharacteristics};
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations requiring arguments already in [0, characteristic - 1]. Not checked.

  /** @brief Same as @ref add, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element add_reduced(Element e1, Element e2) const {
    // e1 + e2 < 2p < 2^bits. For unsigned integers, min(s, s - p) is s - p if s >= p and s otherwise (s - p wraps
    // around): branchless, and compilers turn it into a conditional move.
    const Element s = e1 + e2;
    return std::min<Element>(s, s - characteristic_);
  }

  /** @brief Same as @ref add_inplace, but requires `e1, e2 < characteristic`. */
  void add_inplace_reduced(Element& e1, Element e2) const { e1 = add_reduced(e1, e2); }

  /** @brief Same as @ref subtract, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element subtract_reduced(Element e1, Element e2) const {
    const Element r = e1 - e2;  // wraps around if e1 < e2, then r + p is the result (and r + p < r)
    return std::min<Element>(r, r + characteristic_);
  }

  /** @brief Same as @ref subtract_inplace_front, but requires `e1, e2 < characteristic`. */
  void subtract_inplace_front_reduced(Element& e1, Element e2) const { e1 = subtract_reduced(e1, e2); }

  /** @brief Same as @ref multiply, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element multiply_reduced(Element e1, Element e2) const { return reducer_.mul(e1, e2); }

  /** @brief Same as @ref multiply_inplace, but requires `e1, e2 < characteristic`. */
  void multiply_inplace_reduced(Element& e1, Element e2) const { e1 = multiply_reduced(e1, e2); }

  /** @brief Same as @ref multiply_and_add, but requires `e, m, a < characteristic`. */
  [[nodiscard]] Element multiply_and_add_reduced(Element e, Element m, Element a) const {
    return reducer_.mul_add(e, m, a);
  }

  /** @brief Same as @ref multiply_and_add_inplace_front, but requires `e, m, a < characteristic`. */
  void multiply_and_add_inplace_front_reduced(Element& e, Element m, Element a) const {
    e = multiply_and_add_reduced(e, m, a);
  }

  /** @brief Same as @ref multiply_and_add_inplace_back, but requires `e, m, a < characteristic`. */
  void multiply_and_add_inplace_back_reduced(Element e, Element m, Element& a) const {
    a = multiply_and_add_reduced(e, m, a);
  }

  /** @brief Same as @ref get_inverse, but requires `e < characteristic`. */
  [[nodiscard]] Element get_inverse_reduced(Element e) const {
    return inverse_.empty() ? Gudhi::numbers::mod_inverse<Element, std::make_signed_t<Element> >(e, characteristic_)
                            : inverse_[e];
  }

  /** @brief Same as @ref get_partial_inverse, but requires `e < characteristic`. */
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse_reduced(
      Element e, Characteristic productOfCharacteristics) const {
    return {get_inverse_reduced(e), productOfCharacteristics};
  }

  // ----------------------------------------------------------------------------------------------------------
  // Multiplication by a fixed element (Shoup's method).

  /**
   * @brief Pre-computed data to multiply many elements by the same element. See @ref prepare_multiplier.
   */
  struct Prepared_multiplier {
    Element value;    /**< The multiplier, reduced. */
    Element quotient; /**< `floor(value * 2^bits / characteristic)`. */
  };

  /**
   * @brief Precomputes what is needed to multiply many elements by @p w with @ref multiply_prepared. This is
   * typically useful in loops like `column_i += w * column_j`: a multiplication then costs two integer
   * multiplications and one conditional subtraction.
   */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Prepared_multiplier prepare_multiplier(I w) const {
    const Element v = get_value(w);
    return {v, reducer_.shoup_quotient(v)};
  }

  /** @brief Returns `(x * w.value) % characteristic`. @p x can be any integer, reduced or not. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Element multiply_prepared(I x, const Prepared_multiplier& w) const {
    if constexpr (std::is_unsigned_v<I> && sizeof(I) <= sizeof(Element)) {
      return reducer_.mul_shoup(static_cast<Element>(x), w.value, w.quotient);  // works for any element
    } else {
      return reducer_.mul_shoup(get_value(x), w.value, w.quotient);
    }
  }

  /** @brief Returns `(x * w.value + a) % characteristic`. @p x and @p a can be any integers. */
  template <typename I1, typename I2, class = isIntegers<I1, I2> >
  [[nodiscard]] Element multiply_and_add_prepared(I1 x, const Prepared_multiplier& w, I2 a) const {
    return add_reduced(multiply_prepared(x, w), get_value(a));
  }

  // ----------------------------------------------------------------------------------------------------------

  /** @brief Returns the additive identity of the field, that is 0. */
  static constexpr Element get_additive_identity() { return 0; }

  /** @brief Returns the multiplicative identity of the field, that is 1. */
  static constexpr Element get_multiplicative_identity() { return 1; }

  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   *
   * @param productOfCharacteristics Some value.
   * @return 1.
   */
  static constexpr Element get_partial_multiplicative_identity(
      [[maybe_unused]] Characteristic productOfCharacteristics) {
    return 1;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Zp_field_operators& f1, Zp_field_operators& f2) noexcept {
    std::swap(f1.characteristic_, f2.characteristic_);
    std::swap(f1.reducer_, f2.reducer_);
    f1.inverse_.swap(f2.inverse_);
  }

 private:
  using Reducer = Gudhi::numbers::Modular_reducer<Element>;
  static constexpr unsigned int kBits = sizeof(Element) * 8;
  // Above this characteristic, inverses are computed on demand instead of being stored.
  static constexpr Element kMaxTableCharacteristic = 1U << 16;

  Characteristic characteristic_; /**< Current characteristic of the field. */
  Reducer reducer_;               /**< Fast reductions modulo the characteristic. */
  std::vector<Element> inverse_;  /**< All inverse elements, empty if the characteristic is too large. */

  template <typename U>
  [[nodiscard]] Element _reduce_unsigned(U u) const {
    static_assert(sizeof(U) <= 8, "Integers of more than 64 bits are not supported.");
    return GUDHI_ZP_LIKELY(u < characteristic_) ? static_cast<Element>(u) : reducer_.reduce(u);
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#undef GUDHI_ZP_LIKELY

#endif  // GUDHI_MATRIX_FIELD_ZP_OPERATOR_H_