/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_
#define CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_

namespace Gudhi {

namespace Persistence_representations {

/** \brief The concept Topological_data_with_scalar_product describes the requirements
  * for a type to implement a container that allows computations of scalar products.
  */
class Topological_data_with_scalar_product {
 public:
  double compute_scalar_product(const Topological_data_with_scalar_product& second);
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_
