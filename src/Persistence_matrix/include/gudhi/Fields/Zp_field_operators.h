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
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

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

namespace detail {
// Unsigned integer type with twice the width of the field element type. Used for the intermediate results.
template <std::size_t N>
struct Zp_wide;
template <>
struct Zp_wide<4> {
  using type = std::uint64_t;
};
#ifdef __SIZEOF_INT128__
__extension__ typedef unsigned __int128 Zp_uint128;
template <>
struct Zp_wide<8> {
  using type = Zp_uint128;
};
#endif
}  // namespace detail

/**
 * @class Zp_field_operators Zp_field_operators.h gudhi/Fields/Zp_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 *
 * Design notes:
 * - Elements are always stored fully reduced, that is in \f$ [0, p-1] \f$. Methods without suffix accept any value
 *   and reduce their arguments first. Methods with the suffix `_reduced` require all their arguments to already be
 *   in \f$ [0, p-1] \f$ (not checked) and skip the reduction: use them in hot loops where this is an invariant.
 * - All intermediate results are computed in a type of twice the width of @ref Element, so that nothing overflows.
 * - No hardware division is used for 32 bits elements: reductions use a precomputed Barrett constant.
 *   (64 bits elements use `unsigned __int128` and the compiler's division.)
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
#ifndef __SIZEOF_INT128__
  static_assert(sizeof(Unsigned_integer_type) == 4, "64 bits elements need compiler support for __int128.");
#endif

 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;        /**< Type for the field characteristic. */

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
  Zp_field_operators(Characteristic characteristic = nullCharacteristic)
      : characteristic_(nullCharacteristic), barrett_(0) {
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
    if (!_is_prime(characteristic)) throw std::invalid_argument("Characteristic must be a prime number.");

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
    barrett_ = sizeof(Element) == 4 ? UINT64_MAX / characteristic : 0;
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
  [[nodiscard]] Element add(Element e1, Element e2) const { return add_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Stores in the first element `(e1 + e2) % characteristic`, such that the result is positive. */
  void add_inplace(Element& e1, Element e2) const { e1 = add_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Returns `(e1 - e2) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element subtract(Element e1, Element e2) const { return subtract_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Stores in the first element `(e1 - e2) % characteristic`, such that the result is positive. */
  void subtract_inplace_front(Element& e1, Element e2) const { e1 = subtract_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Stores in the second element `(e1 - e2) % characteristic`, such that the result is positive. */
  void subtract_inplace_back(Element e1, Element& e2) const { e2 = subtract_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Returns `(e1 * e2) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element multiply(Element e1, Element e2) const { return multiply_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Stores in the first element `(e1 * e2) % characteristic`, such that the result is positive. */
  void multiply_inplace(Element& e1, Element e2) const { e1 = multiply_reduced(_reduce(e1), _reduce(e2)); }

  /** @brief Returns `(e * m + a) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element multiply_and_add(Element e, Element m, Element a) const {
    return multiply_and_add_reduced(_reduce(e), _reduce(m), _reduce(a));
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p e. */
  void multiply_and_add_inplace_front(Element& e, Element m, Element a) const {
    e = multiply_and_add_reduced(_reduce(e), _reduce(m), _reduce(a));
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p a. */
  void multiply_and_add_inplace_back(Element e, Element m, Element& a) const {
    a = multiply_and_add_reduced(_reduce(e), _reduce(m), _reduce(a));
  }

  /** @brief Returns `((e + a) * m) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element add_and_multiply(Element e, Element a, Element m) const {
    return multiply_reduced(add_reduced(_reduce(e), _reduce(a)), _reduce(m));
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p e. */
  void add_and_multiply_inplace_front(Element& e, Element a, Element m) const {
    e = multiply_reduced(add_reduced(_reduce(e), _reduce(a)), _reduce(m));
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p m. */
  void add_and_multiply_inplace_back(Element e, Element a, Element& m) const {
    m = multiply_reduced(add_reduced(_reduce(e), _reduce(a)), _reduce(m));
  }

  /** @brief Returns true if `e1 % characteristic == e2 % characteristic`. */
  [[nodiscard]] bool are_equal(Element e1, Element e2) const { return _reduce(e1) == _reduce(e2); }

  /** @brief Returns the opposite of `e % characteristic` in the field. */
  [[nodiscard]] Element get_opposite(Element e) const {
    const Element v = _reduce(e);
    return v == 0 ? Element(0) : characteristic_ - v;
  }

  /** @brief Replaces @p e by its opposite in the field. */
  void get_opposite_inplace(Element& e) const { e = get_opposite(e); }

  /** @brief Returns the inverse of `e % characteristic` in the field. The inverse of 0 is 0 by convention. */
  [[nodiscard]] Element get_inverse(Element e) const { return get_inverse_reduced(_reduce(e)); }

  /** @brief For interface purposes with multi-fields. Returns the inverse together with the second argument. */
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse(Element e,
                                                                       Characteristic productOfCharacteristics) const {
    return {get_inverse(e), productOfCharacteristics};
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations requiring arguments already in [0, characteristic - 1]. Not checked.

  /** @brief Same as @ref add, but requires `e1, e2 < characteristic`. Branchless and overflow free. */
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
  [[nodiscard]] Element multiply_reduced(Element e1, Element e2) const { return _mod_wide(static_cast<Wide>(e1) * e2); }

  /** @brief Same as @ref multiply_inplace, but requires `e1, e2 < characteristic`. */
  void multiply_inplace_reduced(Element& e1, Element e2) const { e1 = multiply_reduced(e1, e2); }

  /** @brief Same as @ref multiply_and_add, but requires `e, m, a < characteristic`. */
  [[nodiscard]] Element multiply_and_add_reduced(Element e, Element m, Element a) const {
    return _mod_wide(static_cast<Wide>(e) * m + a);
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
    return inverse_.empty() ? _inverse_euclid(e, characteristic_) : inverse_[e];
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
  [[nodiscard]] Prepared_multiplier prepare_multiplier(Element w) const {
    w = _reduce(w);
    return {w, static_cast<Element>((static_cast<Wide>(w) << kBits) / characteristic_)};
  }

  /** @brief Returns `(x * w.value) % characteristic`. @p x can be any element, reduced or not. */
  [[nodiscard]] Element multiply_prepared(Element x, const Prepared_multiplier& w) const {
    const Element q = static_cast<Element>((static_cast<Wide>(x) * w.quotient) >> kBits);
    const Element r = static_cast<Element>(x * w.value - q * characteristic_);  // exact value is in [0, 2p)
    return std::min<Element>(r, r - characteristic_);
  }

  /** @brief Returns `(x * w.value + a) % characteristic`. @p x can be any element, @p a has to be reduced. */
  [[nodiscard]] Element multiply_and_add_prepared(Element x, const Prepared_multiplier& w, Element a) const {
    return add_reduced(multiply_prepared(x, w), a);
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
    std::swap(f1.barrett_, f2.barrett_);
    f1.inverse_.swap(f2.inverse_);
  }

 private:
  using Wide = typename detail::Zp_wide<sizeof(Unsigned_integer_type)>::type;
  static constexpr unsigned int kBits = sizeof(Element) * 8;
  // Above this characteristic, inverses are computed on demand instead of being stored.
  static constexpr Element kMaxTableCharacteristic = 1u << 16;

  Characteristic characteristic_; /**< Current characteristic of the field. */
  std::uint64_t barrett_;         /**< floor((2^64 - 1) / characteristic), only used for 32 bits elements. */
  std::vector<Element> inverse_;  /**< All inverse elements, empty if the characteristic is too large. */

  // Barrett reduction of any 64 bits value, 32 bits elements only. At most one correction is needed.
  [[nodiscard]] Element _barrett32(std::uint64_t x) const {
#ifdef __SIZEOF_INT128__
    const std::uint64_t q = static_cast<std::uint64_t>((static_cast<detail::Zp_uint128>(x) * barrett_) >> 64);
    const std::uint64_t r = x - q * characteristic_;  // in [0, 2p)
    return static_cast<Element>(std::min<std::uint64_t>(r, r - characteristic_));
#else
    return static_cast<Element>(x % characteristic_);
#endif
  }

  // Reduces a value of twice the element width, in particular any value < characteristic^2.
  [[nodiscard]] Element _mod_wide(Wide x) const {
    if constexpr (sizeof(Element) == 4) {
      return _barrett32(x);
    } else {
      return static_cast<Element>(x % characteristic_);
    }
  }

  [[nodiscard]] Element _reduce(Element e) const {
    if (GUDHI_ZP_LIKELY(e < characteristic_)) return e;
    if constexpr (sizeof(Element) == 4) {
      return _barrett32(e);
    } else {
      return e % characteristic_;
    }
  }

  template <typename U>
  [[nodiscard]] Element _reduce_unsigned(U u) const {
    if constexpr (sizeof(U) <= sizeof(Element)) {
      return _reduce(static_cast<Element>(u));
    } else if constexpr (sizeof(Element) == 4 && sizeof(U) == 8) {
      return GUDHI_ZP_LIKELY(u < characteristic_) ? static_cast<Element>(u) : _barrett32(u);
    } else {
      return static_cast<Element>(u % characteristic_);
    }
  }

  // Inverse by extended Euclid, a in [0, p). Returns 0 for a == 0. Needs p < 2^(bits-1) for the signed coefficients.
  static Element _inverse_euclid(Element a, Element p) {
    using Signed = std::make_signed_t<Element>;
    Signed t0 = 0, t1 = 1;
    Element r0 = p, r1 = a;
    while (r1 != 0) {
      const Element q = r0 / r1;
      const Element r2 = r0 - q * r1;
      r0 = r1;
      r1 = r2;
      const Signed t2 = t0 - static_cast<Signed>(q) * t1;
      t0 = t1;
      t1 = t2;
    }
    return t0 < 0 ? static_cast<Element>(t0 + static_cast<Signed>(p)) : static_cast<Element>(t0);
  }

  static Element _mulmod_slow(Element a, Element b, Element m) {
    return static_cast<Element>((static_cast<Wide>(a) * b) % m);
  }

  static Element _powmod_slow(Element base, Element exp, Element m) {
    Element res = 1;
    base %= m;
    while (exp > 0) {
      if (exp & 1) res = _mulmod_slow(res, base, m);
      base = _mulmod_slow(base, base, m);
      exp >>= 1;
    }
    return res;
  }

  // Deterministic Miller-Rabin, correct for all n < 3.3 * 10^24, so in particular for all 64 bits integers.
  static bool _is_prime(Element n) {
    constexpr Element smallPrimes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (Element q : smallPrimes) {
      if (n == q) return true;
      if (n % q == 0) return false;
    }
    Element d = n - 1;
    unsigned int s = 0;
    while ((d & 1) == 0) {
      d >>= 1;
      ++s;
    }
    for (Element a : smallPrimes) {
      Element x = _powmod_slow(a, d, n);
      if (x == 1 || x == n - 1) continue;
      bool composite = true;
      for (unsigned int r = 1; r < s; ++r) {
        x = _mulmod_slow(x, x, n);
        if (x == n - 1) {
          composite = false;
          break;
        }
      }
      if (composite) return false;
    }
    return true;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#undef GUDHI_ZP_LIKELY

#endif  // GUDHI_MATRIX_FIELD_ZP_OPERATOR_H_