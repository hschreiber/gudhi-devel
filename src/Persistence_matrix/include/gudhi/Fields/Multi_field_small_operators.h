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
#include <stdexcept>
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
 * 32 bits. Otherwise, use @ref Multi_field_operators.
 *
 * Methods with the suffix `_reduced` require all their arguments to already be in \f$ [0, q-1] \f$ (not checked)
 * and skip the reduction of the arguments.
 */
class Multi_field_operators_with_small_characteristics {
 public:
  using Element = std::uint32_t;                             /**< Type for the elements in the field. */
  using Characteristic = Element;                            /**< Type for the field characteristic. */
  inline static const Characteristic nullCharacteristic = 0; /**< Value of a non initialized characteristic. */

  /**
   * @brief Default constructor. The characteristics have to be set with @ref set_characteristic before use.
   */
  Multi_field_operators_with_small_characteristics() : productOfAllCharacteristics_(nullCharacteristic), barrett_(0) {}

  /**
   * @brief Constructor setting the characteristics.
   *
   * @param minCharacteristic Lower bound of the interval of primes.
   * @param maxCharacteristic Upper bound of the interval of primes (included).
   */
  Multi_field_operators_with_small_characteristics(int minCharacteristic, int maxCharacteristic)
      : productOfAllCharacteristics_(nullCharacteristic), barrett_(0) {
    set_characteristic(minCharacteristic, maxCharacteristic);
  }

  /**
   * @brief Sets the characteristics to all prime numbers in the given interval.
   *
   * @throw std::invalid_argument If the interval is not valid, contains no prime number or if the product of
   * the primes does not fit into 32 bits. If the exception is thrown, the field is left unchanged.
   */
  void set_characteristic(int minimum, int maximum) {
    if (maximum < 2) throw std::invalid_argument("Characteristic must be strictly positive");
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (minimum == maximum && !Gudhi::is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");
    minimum = std::max(minimum, 2);

    std::vector<Element> primes;
    std::uint64_t product = 1;
    for (long long int i = minimum; i <= maximum; ++i) {
      if (Gudhi::is_prime(static_cast<int>(i))) {
        product *= static_cast<std::uint64_t>(i);
        if (product > UINT32_MAX)
          throw std::invalid_argument(
              "The product of the characteristics in the given interval does not fit into "
              "32 bits. Use Multi_field_operators instead.");
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
      const Element R = P / p;
      partials[i] = static_cast<Element>(static_cast<std::uint64_t>(R) * _inverse_euclid(R % p, p));
    }

    // Identities of all sub-products, in an open addressing hash table (at most 2^9 entries for 32 bits).
    // identity(Q) = sum of partials[i] over all primes[i] dividing Q.
    const std::size_t nbSubsets = std::size_t(1) << primes.size();
    unsigned int bits = 1;
    while ((std::size_t(1) << bits) < 2 * nbSubsets) ++bits;
    std::vector<std::pair<Element, Element> > table(std::size_t(1) << bits, {0, 0});  // key 0 = empty
    for (std::size_t mask = 0; mask < nbSubsets; ++mask) {                            // mask 0: Q = 1, identity 0
      std::uint64_t q = 1;
      std::uint64_t id = 0;
      for (std::size_t i = 0; i < primes.size(); ++i) {
        if (mask & (std::size_t(1) << i)) {
          q *= primes[i];
          id += partials[i];
        }
      }
      std::size_t h = _hash(static_cast<Element>(q), bits);
      while (table[h].first != 0) h = (h + 1) & (table.size() - 1);
      table[h] = {static_cast<Element>(q), static_cast<Element>(id % P)};
    }

    primes_.swap(primes);
    productOfAllCharacteristics_ = P;
    barrett_ = UINT64_MAX / P;
    partials_.swap(partials);
    identities_.swap(table);
    identitiesBits_ = bits;
  }

  /** @brief Returns the product of all characteristics. */
  [[nodiscard]] const Characteristic& get_characteristic() const { return productOfAllCharacteristics_; }

  /** @brief Returns `e % characteristic`. */
  [[nodiscard]] Element get_value(Element e) const { return _reduce(e); }

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
    return v == 0 ? v : productOfAllCharacteristics_ - v;
  }

  /** @brief Replaces @p e by its opposite in the field. */
  void get_opposite_inplace(Element& e) const { e = get_opposite(e); }

  /** @brief Returns the partial inverse of @p e with respect to all characteristics, see @ref get_partial_inverse. */
  [[nodiscard]] Element get_inverse(Element e) const {
    return get_partial_inverse(e, productOfAllCharacteristics_).first;
  }

  /**
   * @brief Returns the partial inverse of @p e with respect to the characteristics dividing the given product
   * and the product of the characteristics for which @p e is invertible. The partial inverse is the element which,
   * for each characteristic p dividing the second returned value, is the inverse of @p e modulo p, and is 0 modulo
   * the other characteristics. If @p e is not invertible for any of them, returns `{0, 1}`.
   *
   * @param e Element to get the partial inverse from.
   * @param productOfCharacteristics Product of a subset of the characteristics. 0 is interpreted as the product
   * of all characteristics.
   */
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse(Element e,
                                                                       Characteristic productOfCharacteristics) const {
    GUDHI_CHECK(productOfCharacteristics <= productOfAllCharacteristics_,
                std::invalid_argument(
                    "The given product is not the product of a subset of the current Multi-field characteristics."));
    const Characteristic Q =
        productOfCharacteristics == nullCharacteristic ? productOfAllCharacteristics_ : productOfCharacteristics;
    // One extended Euclid gives both g = gcd(e, Q) and, when g == 1, the inverse of e modulo Q.
    const Element a = Q == productOfAllCharacteristics_ ? _reduce(e) : e % Q;
    Element g;
    Element inv = _extended_euclid(a, Q, g);
    if (g == Q) return {0, get_multiplicative_identity()};  // a == 0: partial inverse is 0
    Characteristic QT = Q;
    if (g != 1) {
      QT = Q / g;
      inv = _extended_euclid(a % QT, QT, g);  // a is invertible modulo QT, as Q is squarefree
    }
    if (QT == productOfAllCharacteristics_) return {inv, QT};
    return {multiply_reduced(_identity_of_subproduct(QT), inv), QT};  // inv < QT <= P
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations requiring arguments already in [0, characteristic - 1]. Not checked.

  /** @brief Same as @ref add, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element add_reduced(Element e1, Element e2) const {
    // computed on 64 bits as P can be close to 2^32; min(s, s - P) is branchless (see Zp_field_operators)
    const std::uint64_t s = static_cast<std::uint64_t>(e1) + e2;
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
    return _barrett(static_cast<std::uint64_t>(e1) * e2);
  }

  /** @brief Same as @ref multiply_inplace, but requires `e1, e2 < characteristic`. */
  void multiply_inplace_reduced(Element& e1, Element e2) const { e1 = multiply_reduced(e1, e2); }

  /** @brief Same as @ref multiply_and_add, but requires `e, m, a < characteristic`. */
  [[nodiscard]] Element multiply_and_add_reduced(Element e, Element m, Element a) const {
    return _barrett((static_cast<std::uint64_t>(e) * m) + a);
  }

  /** @brief Same as @ref multiply_and_add_inplace_front, but requires `e, m, a < characteristic`. */
  void multiply_and_add_inplace_front_reduced(Element& e, Element m, Element a) const {
    e = multiply_and_add_reduced(e, m, a);
  }

  // ----------------------------------------------------------------------------------------------------------

  /** @brief Returns the additive identity of the field, that is 0. */
  static constexpr Element get_additive_identity() { return 0; }

  /** @brief Returns the multiplicative identity of the field, that is 1. */
  static constexpr Element get_multiplicative_identity() { return 1; }

  /**
   * @brief Returns the element which is 1 modulo the characteristics dividing the given product and 0 modulo the
   * others. O(1): the values are tabulated.
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
    std::swap(f1.barrett_, f2.barrett_);
    f1.partials_.swap(f2.partials_);
    f1.identities_.swap(f2.identities_);
    std::swap(f1.identitiesBits_, f2.identitiesBits_);
  }

 private:
  std::vector<Element> primes_;                          /**< All characteristics. */
  Characteristic productOfAllCharacteristics_;           /**< Product of all characteristics. */
  std::uint64_t barrett_;                                /**< floor((2^64 - 1) / productOfAllCharacteristics_). */
  std::vector<Element> partials_;                        /**< Identity of each single characteristic. */
  std::vector<std::pair<Element, Element> > identities_; /**< Hash table: sub-product -> partial identity. */
  unsigned int identitiesBits_ = 0;

  static std::size_t _hash(Element q, unsigned int bits) {
    return static_cast<std::size_t>((static_cast<std::uint32_t>(q) * 2654435769u) >> (32 - bits));
  }

  [[nodiscard]] Element _identity_of_subproduct(Element q) const {
    std::size_t h = _hash(q, identitiesBits_);
    while (identities_[h].first != q) {
      // an empty slot (key 0) means that q is not a product of a subset of the characteristics
      GUDHI_CHECK(identities_[h].first != 0, std::invalid_argument("Not a product of a subset of characteristics."));
      if (identities_[h].first == 0) return 0;
      h = (h + 1) & (identities_.size() - 1);
    }
    return identities_[h].second;
  }

  // Barrett reduction of any 64 bits value. At most one correction is needed.
  [[nodiscard]] Element _barrett(std::uint64_t x) const {
#ifdef __SIZEOF_INT128__
    __extension__ using uint128 = unsigned __int128;
    const auto q = static_cast<std::uint64_t>((static_cast<uint128>(x) * barrett_) >> 64);
    const std::uint64_t r = x - (q * productOfAllCharacteristics_);  // in [0, 2P)
    return static_cast<Element>(std::min<std::uint64_t>(r, r - productOfAllCharacteristics_));
#else
    return static_cast<Element>(x % productOfAllCharacteristics_);
#endif
  }

  [[nodiscard]] Element _reduce(Element e) const {
    return GUDHI_MFS_LIKELY(e < productOfAllCharacteristics_) ? e : _barrett(e);
  }

  // Extended Euclid on (a, m), with a < m, m >= 2: sets g = gcd(a, m) (g = m if a == 0) and returns t in [0, m) such
  // that t * a = g modulo m. In particular, t is the inverse of a modulo m if g == 1 (and is 0 if a == 0).
  static Element _extended_euclid(Element a, Element m, Element& g) {
    std::int64_t t0 = 0, t1 = 1;
    Element r0 = m, r1 = a;
    while (r1 != 0) {
      const Element q = r0 / r1;
      const Element r2 = r0 - (q * r1);
      r0 = r1;
      r1 = r2;
      const std::int64_t t2 = t0 - (static_cast<std::int64_t>(q) * t1);
      t0 = t1;
      t1 = t2;
    }
    g = r0;
    return static_cast<Element>(t0 < 0 ? t0 + static_cast<std::int64_t>(m) : t0);
  }

  static Element _inverse_euclid(Element a, Element m) {
    Element g;
    return _extended_euclid(a, m, g);
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#undef GUDHI_MFS_LIKELY

#endif  // GUDHI_MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_