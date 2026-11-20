/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2026/11 Hannah Schreiber: Factorization with Multi_field_operators
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_
#define PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_

#include <utility>

#include <gudhi/Fields/Multi_field_operators.h>

namespace Gudhi {

namespace persistent_cohomology {

/** \brief Structure representing coefficients in a set of finite fields simultaneously
 * using the chinese remainder theorem.
 *
 * \implements CoefficientField
 * \ingroup persistent_cohomology

 * Details on the algorithms may be found in \cite boissonnat:hal-00922572
 */
class Multi_field {
 public:
  using Operator = persistence_fields::Multi_field_operators;
  using Element = Operator::Element;

  Multi_field() = default;

  /* Initialize the multi-field. The generation of prime numbers might fail with
   * a very small probability.*/
  void init(int min_prime, int max_prime) { op_.set_characteristic(min_prime, max_prime); }

  /** \brief Returns the additive identity \f$0_{\Bbbk}\f$ of the field.*/
  static const Element& additive_identity() { return Operator::get_additive_identity(); }

  /** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
  // 1 everywhere
  static const Element& multiplicative_identity() { return Operator::get_multiplicative_identity(); }

  Element multiplicative_identity(const Element& Q) { return op_.get_partial_multiplicative_identity(Q); }

  /** Returns y * w */
  Element times(const Element& y, const Element& w) { return op_.multiply(y, w); }

  Element plus_equal(const Element& x, const Element& y) { return op_.add(x, y); }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  [[nodiscard]] const Element& characteristic() const { return op_.get_characteristic(); }

  /** Returns the inverse in the field. Modifies P. ??? */
  std::pair<Element, Element> inverse(const Element& x, const Element& QS) { return op_.get_partial_inverse(x, QS); }

  /** Returns -x * y.*/
  Element times_minus(const Element& x, const Element& y) { return op_.multiply(op_.get_opposite(x), y); }

  /** Set x <- x + w * y*/
  Element plus_times_equal(const Element& x, const Element& y, const Element& w) {
    return op_.multiply_and_add(w, y, x);
  }

 private:
  Operator op_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_
