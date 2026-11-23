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
 * @file Multi_field_operators.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_operators class.
 */

#ifndef GUDHI_MATRIX_FIELD_MULTI_OPERATORS_H_
#define GUDHI_MATRIX_FIELD_MULTI_OPERATORS_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <stdexcept>

#include <gmpxx.h>

#include <gudhi/Debug_utils.h>
#include <gudhi/arithmetic.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_operators Multi_field_operators.h gudhi/Fields/Multi_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the multi-field \f$ \mathbb{Z}_{q} \f$ with
 * \f$ q = p_1 \cdots p_k \f$ for all prime numbers \f$ p_i \f$ in a given interval, with arbitrary large \f$ q \f$
 * (using GMP). If \f$ q \f$ fits into 32 bits, @ref Multi_field_operators_with_small_characteristics is much faster.
 *
 * Elements are always stored reduced, that is in \f$ [0, q-1] \f$. Methods without suffix accept any integer, also
 * negative ones. Methods with the suffix `_reduced` require all their arguments to be in \f$ [0, q-1] \f$ (not
 * checked) and avoid some comparisons.
 */
class Multi_field_operators {
 public:
  using Element = mpz_class;                                 /**< Type for the elements in the field. */
  using Characteristic = Element;                            /**< Type for the field characteristic. */
  inline static const Characteristic nullCharacteristic = 0; /**< Value of a non initialized characteristic. */

  /**
   * @brief Default constructor. The characteristics have to be set with @ref set_characteristic before use.
   */
  Multi_field_operators() : productOfAllCharacteristics_(nullCharacteristic), doubleProduct_(nullCharacteristic) {}

  /**
   * @brief Constructor setting the characteristics.
   *
   * @param minCharacteristic Lower bound of the interval of primes.
   * @param maxCharacteristic Upper bound of the interval of primes (included).
   */
  Multi_field_operators(int minCharacteristic, int maxCharacteristic)
      : productOfAllCharacteristics_(nullCharacteristic), doubleProduct_(nullCharacteristic) {
    set_characteristic(minCharacteristic, maxCharacteristic);
  }

  /**
   * @brief Sets the characteristics to all prime numbers in the given interval.
   *
   * @throw std::invalid_argument If the interval is not valid or contains no prime number.
   * If the exception is thrown, the field is left unchanged.
   */
  void set_characteristic(int minimum, int maximum) {
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (maximum < 2) throw std::invalid_argument("The given interval does not contain a positive prime number.");
    if (minimum == maximum && !Gudhi::is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");
    minimum = std::max(minimum, 2);

    std::vector<unsigned int> primes;
    mpz_class tmp(static_cast<unsigned long>(minimum));
    if (mpz_probab_prime_p(tmp.get_mpz_t(), 25) == 0) mpz_nextprime(tmp.get_mpz_t(), tmp.get_mpz_t());
    while (tmp <= maximum) {
      primes.push_back(static_cast<unsigned int>(tmp.get_ui()));
      mpz_nextprime(tmp.get_mpz_t(), tmp.get_mpz_t());
    }
    if (primes.empty()) throw std::invalid_argument("The given interval does not contain a prime number.");

    mpz_class product = 1;
    for (const unsigned int p : primes) product *= p;

    primes_.swap(primes);
    productOfAllCharacteristics_ = product;
    doubleProduct_ = 2 * product;
  }

  /** @brief Returns the product of all characteristics. */
  [[nodiscard]] const Characteristic& get_characteristic() const { return productOfAllCharacteristics_; }

  /** @brief Returns `e % characteristic`, such that the result is positive. */
  [[nodiscard]] Element get_value(Element e) const {
    get_value_inplace(e);
    return e;
  }

  /**
   * @brief Replaces @p e by `e % characteristic`, such that the result is positive.
   */
  void get_value_inplace(Element& e) const {
    if (mpz_sgn(_z(e)) < 0) {
      mpz_add(_z(e), _z(e), _z(productOfAllCharacteristics_));
      if (mpz_sgn(_z(e)) >= 0) return;
    } else {
      if (mpz_cmp(_z(e), _z(productOfAllCharacteristics_)) < 0) return;
      if (mpz_cmp(_z(e), _z(doubleProduct_)) < 0) {
        mpz_sub(_z(e), _z(e), _z(productOfAllCharacteristics_));
        return;
      }
    }
    mpz_mod(_z(e), _z(e), _z(productOfAllCharacteristics_));
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations accepting arbitrary arguments.

  /** @brief Returns `(e1 + e2) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element add(const Element& e1, const Element& e2) const {
    Element r;
    mpz_add(_z(r), _z(e1), _z(e2));
    get_value_inplace(r);
    return r;
  }

  /** @brief Stores in the first element `(e1 + e2) % characteristic`, such that the result is positive. */
  void add_inplace(Element& e1, const Element& e2) const {
    mpz_add(_z(e1), _z(e1), _z(e2));
    get_value_inplace(e1);
  }

  /** @brief Returns `(e1 - e2) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element subtract(const Element& e1, const Element& e2) const {
    Element r;
    mpz_sub(_z(r), _z(e1), _z(e2));
    get_value_inplace(r);
    return r;
  }

  /** @brief Stores in the first element `(e1 - e2) % characteristic`, such that the result is positive. */
  void subtract_inplace_front(Element& e1, const Element& e2) const {
    mpz_sub(_z(e1), _z(e1), _z(e2));
    get_value_inplace(e1);
  }

  /** @brief Stores in the second element `(e1 - e2) % characteristic`, such that the result is positive. */
  void subtract_inplace_back(const Element& e1, Element& e2) const {
    mpz_sub(_z(e2), _z(e1), _z(e2));
    get_value_inplace(e2);
  }

  /** @brief Returns `(e1 * e2) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element multiply(const Element& e1, const Element& e2) const {
    Element r;
    mpz_mul(_z(r), _z(e1), _z(e2));
    get_value_inplace(r);
    return r;
  }

  /** @brief Stores in the first element `(e1 * e2) % characteristic`, such that the result is positive. */
  void multiply_inplace(Element& e1, const Element& e2) const {
    mpz_mul(_z(e1), _z(e1), _z(e2));
    get_value_inplace(e1);
  }

  /** @brief Returns `(e * m + a) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element multiply_and_add(const Element& e, const Element& m, const Element& a) const {
    Element r;
    mpz_mul(_z(r), _z(e), _z(m));
    mpz_add(_z(r), _z(r), _z(a));
    get_value_inplace(r);
    return r;
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p e. */
  void multiply_and_add_inplace_front(Element& e, const Element& m, const Element& a) const {
    mpz_mul(_z(e), _z(e), _z(m));
    mpz_add(_z(e), _z(e), _z(a));
    get_value_inplace(e);
  }

  /** @brief Stores `(e * m + a) % characteristic` in @p a. */
  void multiply_and_add_inplace_back(const Element& e, const Element& m, Element& a) const {
    mpz_addmul(_z(a), _z(e), _z(m));
    get_value_inplace(a);
  }

  /** @brief Returns `((e + a) * m) % characteristic`, such that the result is positive. */
  [[nodiscard]] Element add_and_multiply(const Element& e, const Element& a, const Element& m) const {
    Element r;
    mpz_add(_z(r), _z(e), _z(a));
    mpz_mul(_z(r), _z(r), _z(m));
    get_value_inplace(r);
    return r;
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p e. */
  void add_and_multiply_inplace_front(Element& e, const Element& a, const Element& m) const {
    mpz_add(_z(e), _z(e), _z(a));
    mpz_mul(_z(e), _z(e), _z(m));
    get_value_inplace(e);
  }

  /** @brief Stores `((e + a) * m) % characteristic` in @p m. */
  void add_and_multiply_inplace_back(const Element& e, const Element& a, Element& m) const {
    Element t;
    mpz_add(_z(t), _z(e), _z(a));
    mpz_mul(_z(m), _z(m), _z(t));
    get_value_inplace(m);
  }

  /** @brief Returns true if `e1 % characteristic == e2 % characteristic`. */
  [[nodiscard]] bool are_equal(const Element& e1, const Element& e2) const {
    return mpz_congruent_p(_z(e1), _z(e2), _z(productOfAllCharacteristics_)) != 0;
  }

  /** @brief Returns the opposite of `e % characteristic` in the field. */
  [[nodiscard]] Element get_opposite(const Element& e) const {
    Element v(e);
    get_opposite_inplace(v);
    return v;
  }

  /** @brief Replaces @p e by its opposite in the field. */
  void get_opposite_inplace(Element& e) const {
    get_value_inplace(e);
    if (mpz_sgn(_z(e)) == 0) return;
    mpz_sub(_z(e), _z(productOfAllCharacteristics_), _z(e));
  }

  /** @brief Returns the partial inverse of @p e with respect to all characteristics, see @ref get_partial_inverse. */
  [[nodiscard]] Element get_inverse(const Element& e) const {
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
  [[nodiscard]] std::pair<Element, Characteristic> get_partial_inverse(
      const Element& e, const Characteristic& productOfCharacteristics) const {
    GUDHI_CHECK(productOfCharacteristics >= 0 && productOfCharacteristics <= productOfAllCharacteristics_,
                std::invalid_argument(
                    "The given product is not the product of a subset of the current Multi-field characteristics."));
    const Characteristic& Q =
        mpz_sgn(_z(productOfCharacteristics)) == 0 ? productOfAllCharacteristics_ : productOfCharacteristics;
    Characteristic g;
    mpz_gcd(_z(g), _z(e), _z(Q));
    if (g == Q) return {0, get_multiplicative_identity()};  // partial inverse is 0
    std::pair<Element, Characteristic> res;
    Characteristic& QT = res.second;
    Element& x = res.first;
    mpz_divexact(_z(QT), _z(Q), _z(g));
    if (QT == productOfAllCharacteristics_) {
      mpz_mod(_z(x), _z(e), _z(QT));
      mpz_invert(_z(x), _z(x), _z(QT));
      return res;
    }
    Characteristic R;
    mpz_divexact(_z(R), _z(productOfAllCharacteristics_), _z(QT));
    mpz_mul(_z(x), _z(e), _z(R));
    mpz_mod(_z(x), _z(x), _z(QT));
    mpz_invert(_z(x), _z(x), _z(QT));
    mpz_mul(_z(x), _z(x), _z(R));  // < R * QT = q
    return res;
  }

  // ----------------------------------------------------------------------------------------------------------
  // Operations requiring arguments already in [0, characteristic - 1]. Not checked.

  /** @brief Same as @ref add, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element add_reduced(const Element& e1, const Element& e2) const {
    Element r;
    mpz_add(_z(r), _z(e1), _z(e2));
    if (mpz_cmp(_z(r), _z(productOfAllCharacteristics_)) >= 0) mpz_sub(_z(r), _z(r), _z(productOfAllCharacteristics_));
    return r;
  }

  /** @brief Same as @ref add_inplace, but requires `e1, e2 < characteristic`. */
  void add_inplace_reduced(Element& e1, const Element& e2) const {
    mpz_add(_z(e1), _z(e1), _z(e2));
    if (mpz_cmp(_z(e1), _z(productOfAllCharacteristics_)) >= 0)
      mpz_sub(_z(e1), _z(e1), _z(productOfAllCharacteristics_));
  }

  /** @brief Same as @ref subtract, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element subtract_reduced(const Element& e1, const Element& e2) const {
    Element r;
    mpz_sub(_z(r), _z(e1), _z(e2));
    if (mpz_sgn(_z(r)) < 0) mpz_add(_z(r), _z(r), _z(productOfAllCharacteristics_));
    return r;
  }

  /** @brief Same as @ref subtract_inplace_front, but requires `e1, e2 < characteristic`. */
  void subtract_inplace_front_reduced(Element& e1, const Element& e2) const {
    mpz_sub(_z(e1), _z(e1), _z(e2));
    if (mpz_sgn(_z(e1)) < 0) mpz_add(_z(e1), _z(e1), _z(productOfAllCharacteristics_));
  }

  /** @brief Same as @ref multiply, but requires `e1, e2 < characteristic`. */
  [[nodiscard]] Element multiply_reduced(const Element& e1, const Element& e2) const {
    Element r;
    mpz_mul(_z(r), _z(e1), _z(e2));
    mpz_tdiv_r(_z(r), _z(r), _z(productOfAllCharacteristics_));
    return r;
  }

  /** @brief Same as @ref multiply_inplace, but requires `e1, e2 < characteristic`. */
  void multiply_inplace_reduced(Element& e1, const Element& e2) const {
    mpz_mul(_z(e1), _z(e1), _z(e2));
    mpz_tdiv_r(_z(e1), _z(e1), _z(productOfAllCharacteristics_));
  }

  /** @brief Same as @ref multiply_and_add, but requires `e, m, a < characteristic`. */
  [[nodiscard]] Element multiply_and_add_reduced(const Element& e, const Element& m, const Element& a) const {
    Element r;
    mpz_mul(_z(r), _z(e), _z(m));
    mpz_add(_z(r), _z(r), _z(a));
    mpz_tdiv_r(_z(r), _z(r), _z(productOfAllCharacteristics_));
    return r;
  }

  /** @brief Same as @ref multiply_and_add_inplace_front, but requires `e, m, a < characteristic`. */
  void multiply_and_add_inplace_front_reduced(Element& e, const Element& m, const Element& a) const {
    mpz_mul(_z(e), _z(e), _z(m));
    mpz_add(_z(e), _z(e), _z(a));
    mpz_tdiv_r(_z(e), _z(e), _z(productOfAllCharacteristics_));
  }

  // ----------------------------------------------------------------------------------------------------------

  /** @brief Returns a const reference to the additive identity of the field, that is 0. */
  static const Element& get_additive_identity() { return additiveID_; }

  /** @brief Returns a const reference to the multiplicative identity of the field, that is 1. */
  static const Element& get_multiplicative_identity() { return multiplicativeID_; }

  /**
   * @brief Returns the element which is 1 modulo the characteristics dividing the given product and 0 modulo the
   * others, that is `R * (R^{-1} mod Q)` with `R = q / Q`.
   *
   * @param productOfCharacteristics Product of a subset of the characteristics. 0 is interpreted as the product
   * of all characteristics.
   */
  [[nodiscard]] Element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics) const {
    GUDHI_CHECK(productOfCharacteristics >= 0 && productOfCharacteristics <= productOfAllCharacteristics_,
                std::invalid_argument(
                    "The given product is not the product of a subset of the current Multi-field characteristics."));
    if (productOfCharacteristics == nullCharacteristic || productOfCharacteristics == productOfAllCharacteristics_) {
      return get_multiplicative_identity();
    }
    if (productOfCharacteristics == 1) return get_additive_identity();
    Element R;
    mpz_divexact(_z(R), _z(productOfAllCharacteristics_), _z(productOfCharacteristics));
    Element x;
    mpz_mod(_z(x), _z(R), _z(productOfCharacteristics));
    mpz_invert(_z(x), _z(x), _z(productOfCharacteristics));
    mpz_mul(_z(x), _z(x), _z(R));
    return x;
  }

  friend void swap(Multi_field_operators& f1, Multi_field_operators& f2) noexcept {
    f1.primes_.swap(f2.primes_);
    std::swap(f1.productOfAllCharacteristics_, f2.productOfAllCharacteristics_);
    std::swap(f1.doubleProduct_, f2.doubleProduct_);
  }

 private:
  std::vector<unsigned int> primes_;           /**< All characteristics. */
  Characteristic productOfAllCharacteristics_; /**< Product of all characteristics. */
  Characteristic doubleProduct_;               /**< Twice the product of all characteristics. */
  inline static const Element multiplicativeID_ = 1;
  inline static const Element additiveID_ = 0;

  static mpz_ptr _z(Element& e) { return e.get_mpz_t(); }
  static mpz_srcptr _z(const Element& e) { return e.get_mpz_t(); }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // GUDHI_MATRIX_FIELD_MULTI_OPERATORS_H_