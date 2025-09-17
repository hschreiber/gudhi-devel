/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2025/09 Hannah Schreiber: Replacing Field_Zp with more generic Zp_field_operators
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
#define PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

#include <gudhi/Fields/Zp_field_operators.h>

namespace Gudhi {

namespace persistent_cohomology {

using Field_Zp = Gudhi::persistence_fields::Zp_field_operators<unsigned int>;

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
