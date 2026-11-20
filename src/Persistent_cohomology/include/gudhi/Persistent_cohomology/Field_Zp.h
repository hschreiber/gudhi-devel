/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2026/11 Hannah Schreiber: Factorization with Zp_field_operators
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
#define PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

#include <cstdint>
#include <utility>

#include <gudhi/Fields/Zp_field_operators.h>

namespace Gudhi {

namespace persistent_cohomology {

/** \brief Structure representing the coefficient field \f$\mathbb{Z}/p\mathbb{Z}\f$
 *
 * \implements CoefficientField
 * \ingroup persistent_cohomology
 */
class Field_Zp {
 public:
  using Element = int;
  using Operator = persistence_fields::Zp_field_operators<std::uint32_t>;

  Field_Zp() = default;

  void init(int charac) { op_.set_characteristic(charac); }

  /** Set x <- x + w * y*/
  Element plus_times_equal(const Element& x, const Element& y, const Element& w) {
    return static_cast<Element>(op_.multiply_and_add(w, y, x));
  }

  /** Returns y * w */
  Element times(Element y, Element w) { return static_cast<Element>(op_.multiply(y, w)); }

  Element plus_equal(Element x, Element y) { return static_cast<Element>(op_.add(x, y)); }

  /** \brief Returns the additive identity \f$0_{\Bbbk}\f$ of the field.*/
  static Element additive_identity() { return Operator::get_additive_identity(); }

  /** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
  static Element multiplicative_identity([[maybe_unused]] Element e = 0) {
    return Operator::get_multiplicative_identity();
  }

  /** Returns the inverse in the field. Modifies P. ??? */
  // <------ return the product of field characteristic for which x is invertible
  std::pair<Element, Element> inverse(Element x, Element P) { return op_.get_partial_inverse(x, P); }

  /** Returns -x * y.*/
  Element times_minus(Element x, Element y) { return static_cast<Element>(op_.multiply(op_.get_opposite(x), y)); }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  [[nodiscard]] int characteristic() const { return static_cast<int>(op_.get_characteristic()); }

 private:
  Operator op_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
