/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, Clément Maria
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_field_small_operators.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_operators_with_small_characteristics class.
 */

#ifndef GUDHI_MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_
#define GUDHI_MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/arithmetic.h>

#if defined(__GNUC__) || defined(__clang__)
#define GUDHI_MFS_LIKELY(x) __builtin_expect(!!(x), 1)
#else
#define GUDHI_MFS_LIKELY(x) (x)
#endif

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_operators_with_small_characteristics Multi_field_small_operators.h \
 * gudhi/Fields/Multi_field_small_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the multi-field \f$ \mathbb{Z}_{q} \f$ with
 * \f$ q = p_1 \cdots p_k \f$ for all prime numbers \f$ p_i \f$ in a given interval, as long as \f$ q \f$ fits into
 * the element type (strictly below \f$ 2^{63} \f$ for 64 bits). Otherwise, use @ref Multi_field_operators.
 *
 * All intermediate results are computed on at least 64 bits, such that only \f$ q \f$ has to fit into the element
 * type. For example, with the default 32 bits type, the primes 2 to 23 are possible (instead of
 * 2 to 13), and with `std::uint64_t` the primes 2 to 47.
 *
 * @tparam Unsigned_integer_type A native unsigned integer type of at most 64 bits (`bool` excluded). Will be used
 * as the field element type.
 *
 * Methods without suffix accept any native integer type (signed or not, up to 64 bits, a different type for each
 * argument) and any value, also negative ones, and reduce their arguments first. They return an @ref Element. The
 * "in place" versions store the result in an @ref Element reference (the value it contains before the call can be any
 * @ref Element, reduced or not). Methods with the suffix `_reduced` only accept @ref Element arguments which are
 * already in \f$ [0, q-1] \f$ (no check) and skip the reduction.
 */
template <typename Unsigned_integer_type = std::uint32_t>
class Multi_field_operators_with_small_characteristics {
  static_assert(std::is_unsigned_v<Unsigned_integer_type> && !std::is_same_v<Unsigned_integer_type, bool> &&
                    sizeof(Unsigned_integer_type) <= 8,
                "Multi_field_operators_with_small_characteristics supports only unsigned integer types of at most "
                "64 bits.");

 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;        /**< Type for the field characteristic. */
  template <typename... Ts>
  using isIntegers = std::enable_if_t<(std::is_integral_v<Ts> && ...)>;

 private:
  using Reducer = Gudhi::numbers::Modular_reducer<Element>;
  /**
   * Type used for the intermediate values: 32 bits for types of at most 32 bits, else the type itself.
   * That is, Red_element will coincide with Element only when Element is std::uint32_t or std::uint64_t.
   */
  using Red_element = typename Reducer::Element;

 public:
  inline static const Characteristic nullCharacteristic = 0; /**< Value of a non initialized characteristic. */

  /**
   * @brief Default constructor. The characteristics have to be set with @ref set_characteristic before use.
   */
  Multi_field_operators_with_small_characteristics() : productOfAllCharacteristics_(nullCharacteristic) {}

  /**
   * @brief Constructor setting the characteristics.
   *
   * @param minCharacteristic Lower bound of the interval of primes.
   * @param maxCharacteristic Upper bound of the interval of primes (included).
   */
  Multi_field_operators_with_small_characteristics(int minCharacteristic, int maxCharacteristic)
      : productOfAllCharacteristics_(nullCharacteristic) {
    set_characteristic(minCharacteristic, maxCharacteristic);
  }

  /**
   * @brief Sets the characteristics to all prime numbers in the given interval.
   *
   * @throw std::invalid_argument If the interval is not valid, contains no prime number or if the product of
   * the primes does not fit into the element type. If the exception is thrown, the field is left unchanged.
   */
  void set_characteristic(int minimum, int maximum) {
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (maximum < 2) throw std::invalid_argument("The given interval does not contain a positive prime number.");
    if (minimum == maximum && !Gudhi::is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");
    minimum = std::max(minimum, 2);

    // The product has to fit in the element type, and has to be below 2^63 for 64 bits (for the signed Euclid
    // coefficients and to be able to add two elements without overflow).
    const std::uint64_t limit = sizeof(Element) == 8 ? static_cast<std::uint64_t>(INT64_MAX)
                                                     : static_cast<std::uint64_t>(std::numeric_limits<Element>::max());
    std::vector<Element> primes;
    std::uint64_t product = 1;
    for (long long int i = minimum; i <= maximum; ++i) {
      if (Gudhi::is_prime(static_cast<int>(i))) {
        if (static_cast<std::uint64_t>(i) > limit / product)  // product * i > limit
          throw std::invalid_argument(
              "The product of the characteristics in the given interval does not fit into "
              "the element type. Use a larger type or Multi_field_operators instead.");
        product *= static_cast<std::uint64_t>(i);
        primes.push_back(static_cast<Element>(i));
      }
    }
    if (primes.empty()) throw std::invalid_argument("The given interval does not contain a prime number.");

    const auto P = static_cast<Characteristic>(product);
    // partial[i] is the element which is 1 modulo primes[i] and 0 modulo all other primes: R * (R^{-1} mod p)
    // with R = P / p.
    std::vector<Element> partials(primes.size());
    for (std::size_t i = 0; i < primes.size(); ++i) {
      const Element p = primes[i];
      const auto R = static_cast<Element>(P / p);
      partials[i] =
          static_cast<Element>(static_cast<std::uint64_t>(R) * Gudhi::numbers::mod_inverse<Element>(R % p, p));
    }

    // Identities of all sub-products, in an open addressing hash table (at most 2^12 entries, and none above, where
    // identities are computed on demand).
    // identity(Q) = sum of partials[i] over all primes[i] dividing Q.
    unsigned int bits = 0;
    std::vector<std::pair<Element, Element> > table;
    if (primes.size() <= kMaxTabulatedCharacteristics) {
      const std::size_t nbSubsets = std::size_t(1) << primes.size();
      bits = 1;
      while ((std::size_t(1) << bits) < 2 * nbSubsets) ++bits;
      table.assign(std::size_t(1) << bits, {0, 0});           // key 0 = empty
      for (std::size_t mask = 0; mask < nbSubsets; ++mask) {  // mask 0: Q = 1, identity 0
        std::uint64_t q = 1;
        std::uint64_t id = 0;
        for (std::size_t i = 0; i < primes.size(); ++i) {
          if ((mask & (std::size_t(1) << i)) != 0U) {
            q *= primes[i];
            const std::uint64_t pt = partials[i];
            id = id >= P - pt ? id - (P - pt) : id + pt;  // modular addition, without overflow
          }
        }
        std::size_t h = _hash(static_cast<Element>(q), bits);
        while (table[h].first != 0) h = (h + 1) & (table.size() - 1);
        table[h] = {static_cast<Element>(q), static_cast<Element>(id)};
      }
    }

    primes_.swap(primes);
    productOfAllCharacteristics_ = P;
    reducer_ = Reducer(P);
    partials_.swap(partials);
    identities_.swap(table);
    identitiesBits_ = bits;
  }

  /** @brief Returns the product of all characteristics. */
  [[nodiscard]] const Characteristic& get_characteristic() const { return productOfAllCharacteristics_; }

  /**
   * @brief Returns the value of an integer in the field: `e % characteristic`, such that the result is positive.
   *
   * @tparam Integer_type Any native integer type, signed or not, of at most 64 bits.
   */
  template <typename Integer_type, class = isIntegers<Integer_type> >
  [[nodiscard]] Element get_value(Integer_type e) const {
    if constexpr (std::is_signed_v<Integer_type>) {
      using U = std::make_unsigned_t<Integer_type>;
      if (e >= 0) return _reduce_unsigned(static_cast<U>(e));
      const Element r = _reduce_unsigned(static_cast<U>(U(0) - static_cast<U>(e)));  // reduces |e|
      return r == 0 ? Element(0) : static_cast<Element>(productOfAllCharacteristics_ - r);
    } else {
      return _reduce_unsigned(e);
    }
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations accepting any integer type and any value.

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
    return v == 0 ? v : static_cast<Element>(productOfAllCharacteristics_ - v);
  }

  /** @brief Replaces @p e by its opposite in the field. */
  void get_opposite_inplace(Element& e) const { e = get_opposite(e); }

  /** @brief Returns the partial inverse of @p e with respect to all characteristics, see @ref get_partial_inverse. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Element get_inverse(I e) const {
    return get_partial_inverse(e, productOfAllCharacteristics_).first;
  }

  /**
   * @brief Returns the partial inverse of @p e with respect to the characteristics dividing the given product
   * and the product of the characteristics for which @p e is invertible. The partial inverse is the element which,
   * for each characteristic p dividing the second returned value, is the inverse of @p e modulo p, and is 0 modulo
   * the other characteristics. If @p e is not invertible for any of them, returns `{0, 1}`.
   *
   * @param e Integer to get the partial inverse from, of any native integer type.
   * @param productOfCharacteristics Product of a subset of the characteristics. 0 is interpreted as the product
   * of all characteristics.
   */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse(I e,
                                                                       Characteristic productOfCharacteristics) const {
    GUDHI_CHECK(productOfCharacteristics <= productOfAllCharacteristics_,
                std::invalid_argument(
                    "The given product is not the product of a subset of the current Multi-field characteristics."));
    const Characteristic Q =
        productOfCharacteristics == nullCharacteristic ? productOfAllCharacteristics_ : productOfCharacteristics;
    // One extended Euclid gives both g = gcd(e, Q) and, when g == 1, the inverse of e modulo Q.
    // (e modulo the product, then modulo Q, is e modulo Q as Q divides the product.)
    const Element v = get_value(e);
    const Element a = Q == productOfAllCharacteristics_ ? v : static_cast<Element>(v % Q);
    Element g;
    Element inv = Gudhi::numbers::extended_euclid<Element>(a, Q, g);
    if (g == Q) return {0, get_multiplicative_identity()};  // a == 0: partial inverse is 0
    Characteristic QT = Q;
    if (g != 1) {
      QT = static_cast<Characteristic>(Q / g);
      // a is invertible modulo QT, as Q is squarefree
      inv = Gudhi::numbers::extended_euclid<Element>(static_cast<Element>(a % QT), QT, g);
    }
    if (QT == productOfAllCharacteristics_) return {inv, QT};
    return {multiply_reduced(_identity_of_subproduct(QT), inv), QT};  // inv < QT <= P
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations requiring arguments already in [0, characteristic - 1]. Not checked.

  /** @brief Same as @ref add, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element add_reduced(Element e1, Element e2) const {
    // computed on 64 bits as P can be close to 2^32; min(s, s - P) is branchless (see Zp_field_operators)
    const std::uint64_t s = static_cast<std::uint64_t>(e1) + e2;  // < 2 * 2^63 for 64 bits elements
    return static_cast<Element>(std::min<std::uint64_t>(s, s - productOfAllCharacteristics_));
  }

  /** @brief Same as @ref add_inplace, but requires `e1, e2 < characteristic`. */
  void add_inplace_reduced(Element& e1, Element e2) const { e1 = add_reduced(e1, e2); }

  /** @brief Same as @ref subtract, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element subtract_reduced(Element e1, Element e2) const {
    const std::uint64_t r = static_cast<std::uint64_t>(e1) + productOfAllCharacteristics_ - e2;  // in (0, 2P)
    return static_cast<Element>(std::min<std::uint64_t>(r, r - productOfAllCharacteristics_));
  }

  /** @brief Same as @ref subtract_inplace_front, but requires `e1, e2 < characteristic`. */
  void subtract_inplace_front_reduced(Element& e1, Element e2) const { e1 = subtract_reduced(e1, e2); }

  /** @brief Same as @ref multiply, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element multiply_reduced(Element e1, Element e2) const {
    return static_cast<Element>(reducer_.mul(e1, e2));
  }
  void multiply_inplace_reduced(Element& e1, Element e2) const { e1 = multiply_reduced(e1, e2); }
  [[nodiscard]] Element multiply_and_add_reduced(Element e, Element m, Element a) const {
    return static_cast<Element>(reducer_.mul_add(e, m, a));
  }
  void multiply_and_add_inplace_front_reduced(Element& e, Element m, Element a) const {
    e = multiply_and_add_reduced(e, m, a);
  }

  // ----------------------------------------------------------------------------------------------------------
  // Multiplication by a fixed element (Shoup's method).

  /** @brief Pre-computed data to multiply many elements by the same element. See @ref prepare_multiplier. */
  struct Prepared_multiplier {
    Element value;        /**< The multiplier, reduced. */
    Red_element quotient; /**< `floor(value * 2^b / characteristic)`, with b = 32 for types of at most 32 bits. */
  };

  /**
   * @brief Pre-computes what is needed to multiply many elements by @p w with @ref multiply_prepared, typically in
   * loops like `column_i += w * column_j`: a multiplication then costs a few integer multiplications and one
   * conditional subtraction, and no reduction by the characteristic.
   */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Prepared_multiplier prepare_multiplier(I w) const {
    const Element v = get_value(w);
    return {v, reducer_.shoup_quotient(v)};
  }

  /** @brief Returns `(x * w.value) % characteristic`. @p x can be any integer, reduced or not. */
  template <typename I, class = isIntegers<I> >
  [[nodiscard]] Element multiply_prepared(I x, const Prepared_multiplier& w) const {
    if constexpr (std::is_unsigned_v<I> && sizeof(I) <= sizeof(Red_element)) {
      return static_cast<Element>(reducer_.mul_shoup(static_cast<Red_element>(x), w.value, w.quotient));  // any x works
    } else {
      return static_cast<Element>(reducer_.mul_shoup(get_value(x), w.value, w.quotient));
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
   * @brief Returns the element which is 1 modulo the characteristics dividing the given product and 0 modulo the
   * others.
   *
   * @param productOfCharacteristics Product of a subset of the characteristics. 0 is interpreted as the product
   * of all characteristics.
   */
  [[nodiscard]] Element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics) const {
    GUDHI_CHECK(productOfCharacteristics <= productOfAllCharacteristics_,
                std::invalid_argument(
                    "The given product is not the product of a subset of the current Multi-field characteristics."));
    if (productOfCharacteristics == nullCharacteristic || productOfCharacteristics == productOfAllCharacteristics_) {
      return get_multiplicative_identity();
    }
    return _identity_of_subproduct(productOfCharacteristics);
  }

  friend void swap(Multi_field_operators_with_small_characteristics& f1,
                   Multi_field_operators_with_small_characteristics& f2) noexcept {
    f1.primes_.swap(f2.primes_);
    std::swap(f1.productOfAllCharacteristics_, f2.productOfAllCharacteristics_);
    std::swap(f1.reducer_, f2.reducer_);
    f1.partials_.swap(f2.partials_);
    f1.identities_.swap(f2.identities_);
    std::swap(f1.identitiesBits_, f2.identitiesBits_);
  }

 private:
  std::vector<Element> primes_;                          /**< All characteristics. */
  Characteristic productOfAllCharacteristics_;           /**< Product of all characteristics. */
  static constexpr std::size_t kMaxTabulatedCharacteristics = 12;
  Reducer reducer_;                                      /**< Fast reductions modulo the product. */
  std::vector<Element> partials_;                        /**< Identity of each single characteristic. */
  std::vector<std::pair<Element, Element> > identities_; /**< Hash table: sub-product -> partial identity. */
  unsigned int identitiesBits_ = 0;

  static std::size_t _hash(Element q, unsigned int bits) {
    if constexpr (sizeof(Element) == 8) {
      return static_cast<std::size_t>((static_cast<std::uint64_t>(q) * 0x9E3779B97F4A7C15ULL) >> (64 - bits));
    } else {
      return static_cast<std::size_t>((static_cast<std::uint32_t>(q) * 2654435769U) >> (32 - bits));
    }
  }

  [[nodiscard]] Element _identity_of_subproduct(Element q) const {
    if (identities_.empty()) {  // too many characteristics for a table: R * (R^{-1} mod q), with R = P / q
      if (q == 1) return 0;
      const auto R = static_cast<Element>(productOfAllCharacteristics_ / q);
      return static_cast<Element>(static_cast<std::uint64_t>(R) *
                                  Gudhi::numbers::mod_inverse<Element>(static_cast<Element>(R % q), q));
    }
    std::size_t h = _hash(q, identitiesBits_);
    while (identities_[h].first != q) {
      // an empty slot (key 0) means that q is not a product of a subset of the characteristics
      GUDHI_CHECK(identities_[h].first != 0, std::invalid_argument("Not a product of a subset of characteristics."));
      if (identities_[h].first == 0) return 0;
      h = (h + 1) & (identities_.size() - 1);
    }
    return identities_[h].second;
  }

  template <typename U>
  [[nodiscard]] Element _reduce_unsigned(U u) const {
    static_assert(sizeof(U) <= 8, "Integers of more than 64 bits are not supported.");
    return GUDHI_MFS_LIKELY(u < productOfAllCharacteristics_) ? static_cast<Element>(u)
                                                              : static_cast<Element>(reducer_.reduce(u));
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#undef GUDHI_MFS_LIKELY

#endif  // MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_