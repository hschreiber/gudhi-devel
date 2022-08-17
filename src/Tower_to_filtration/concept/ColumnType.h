/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  TU Graz (Austria)
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

#ifndef CONCEPT_COLUMN_TYPE_H_
#define CONCEPT_COLUMN_TYPE_H_

/** @file ColumnType.h
 * @brief Contains @ref Gudhi::tower_to_filtration::ColumnType concept.
 */

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Data structure for the columns of the boundary matrix used by @ref Persistence.
 * This concept describes the methodes required by the module.
 * The considered cofficients are in @f$\mathbb{Z}_2@f$ only. Therefore the values of the columns are the indices of non-zero cells in increasing order.
 */
class ColumnType{
public:
    /**
     * Type for cell content. Should correspond to @ref Gudhi::tower_to_filtration::Persistence::index, if used for @ref Gudhi::tower_to_filtration::Persistence.
     * Needs to be comparable (== and <) and order corresponds to order in filtration.
     */
    typename coefficient_type;

    /**
     * @brief Constructor
     * @param dim dimension of the simplex whose boundary is represented by this column.
     */
    ColumnType(int dim);

    /**
     * @brief Insert a cell at the end of the column.
     * @param cell row index of the cell to be inserted.
     */
    void push_back(coefficient_type cell);
    /**
     * @brief Return the pivot of the column, which is the index of the last nonzero value.
     * @return The last nonzero value of the column.
     */
    coefficient_type get_pivot();
    /**
     * @brief Sum this column and another column. The result replaces this column.
     * @param columnToAdd column to add up.
     */
    void add(ColumnType &columnToAdd);
    /**
     * @brief Erase useless cells from the column. See @cite KerberS17 for more details.
     * @param latest private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
     * @param isActivePositive private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
     * @param columns private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
     */
    void clean(std::unordered_map<coefficient_type, coefficient_type> *latest, std::unordered_map<coefficient_type, std::pair<bool, bool> *> *isActivePositive,
	       std::unordered_map<coefficient_type, List_column *> *columns);
};

}
}

#endif  // CONCEPT_COLUMN_TYPE_H_
