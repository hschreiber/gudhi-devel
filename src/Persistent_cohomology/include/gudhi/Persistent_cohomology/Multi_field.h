/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2025/09 Hannah Schreiber: Replacing Multi_field with more generic Multi_field_operators
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_
#define PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_

#include <gudhi/Fields/Multi_field_operators.h>

namespace Gudhi {

namespace persistent_cohomology {

using Multi_field = Gudhi::persistence_fields::Multi_field_operators;

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_MULTI_FIELD_H_
