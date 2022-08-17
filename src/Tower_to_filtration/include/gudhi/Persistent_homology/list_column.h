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

#ifndef LISTCOLUMN_H
#define LISTCOLUMN_H

/** @file list_column.h
 * @brief Contains @ref Gudhi::tower_to_filtration::List_column class.
 */

#include <iostream>
#include <list>
#include <unordered_map>

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Column type which is based on `std::list<coefficient_type>`. Fulfills the requirements of the @ref Gudhi::tower_to_filtration::ColumnType concept.
 * Its coefficient are in @f$\mathbb{Z}_2@f$ only. Therefore the values of the list are the indices of non-zero cells in increasing order.
 */
class List_column
{
public:
    /**
     * Type for cell content. Should correspond to @ref Gudhi::tower_to_filtration::Persistence::index, if used for @ref Gudhi::tower_to_filtration::Persistence.
     * Needs to be comparable (== and <) and order corresponds to order in filtration.
     */
    using coefficient_type = long long;

    List_column(int dim);
    ~List_column();

    void add(List_column &columnToAdd);
    /**
     * @brief Erase a cell from the column.
     * @param pos position of the cell to be deleted.
     */
    void erase(std::list<coefficient_type>::iterator &pos){ column_->erase(pos); }
    /**
     * @brief Return an iterator of the column pointing at the begining.
     * @return An iterator of the column pointing at the begining.
     */
    std::list<coefficient_type>::iterator get_begin_iterator(){ return column_->begin(); }
    /**
     * @brief Return an reverse iterator of the column pointing at the end.
     * @return An reverse iterator of the column pointing at the end.
     */
    std::list<coefficient_type>::reverse_iterator get_reverse_begin_iterator(){ return column_->rbegin(); }
    /**
     * @brief Return an iterator of the column pointing after the end.
     * @return An iterator of the column pointing after the end.
     */
    std::list<coefficient_type>::iterator get_end_iterator(){ return column_->end(); }
    /**
     * @brief Return an reverse iterator of the column pointing before the begining.
     * @return An reverse iterator of the column pointing before the begining.
     */
    std::list<coefficient_type>::reverse_iterator get_reverse_end_iterator(){ return column_->rend(); }
    /**
     * @brief Returns the number of nonzero values in the column.
     * @return The number of nonzero values in the column.
     */
    std::list<coefficient_type>::size_type get_size(){ return column_->size(); }
    /**
     * @brief Returns the stored dimension.
     * @return The dimension.
     */
    int get_dim() const{ return dim_; }
    coefficient_type get_pivot();
    void clean(std::unordered_map<coefficient_type, coefficient_type> *latest, std::unordered_map<coefficient_type, std::pair<bool, bool> *> *isActivePositive,
	       std::unordered_map<coefficient_type, List_column *> *columns);
    void push_back(coefficient_type cell);

private:
    int dim_;				    /**< Dimension of the column. */
    std::list<coefficient_type> *column_;   /**< Data container of the column. */
};

/**
 * @brief Constructor
 * @param dim dimension to be stored as the column dimension.
 */
inline List_column::List_column(int dim) : dim_(dim)
{
    column_ = new std::list<coefficient_type>();
}

/**
 * @brief Destructor
 */
inline List_column::~List_column()
{
    delete column_;
}

/**
 * @brief Replaces the column values by the sum of this column and @p columnToAdd.
 * @param columnToAdd column to sum with.
 */
inline void List_column::add(List_column &columnToAdd)
{
    std::list<coefficient_type>::iterator itToAdd = columnToAdd.get_begin_iterator(), itTarget = column_->begin();
    while (itToAdd != columnToAdd.get_end_iterator() && itTarget != column_->end()){
        if (*itToAdd == *itTarget){
            column_->erase(itTarget++);
            itToAdd++;
        } else if (*itToAdd < *itTarget){
            column_->insert(itTarget, *itToAdd);
            itToAdd++;
        } else {
            itTarget++;
        }
    }
    while (itToAdd != columnToAdd.get_end_iterator()){
        column_->push_back(*itToAdd);
        itToAdd++;
    }
}

/**
 * @brief Returns the pivot of the column, i.e. the index of the last nonzero value.
 * @return The pivot of the column.
 */
inline List_column::coefficient_type List_column::get_pivot(){
    if (column_->empty()) return -1;
    return column_->back();
}

/**
 * @brief Erase useless cells from the column. See @cite KerberS17 for more details.
 * @param latest private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 * @param isActivePositive private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 * @param columns private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 */
inline void List_column::clean(std::unordered_map<coefficient_type, coefficient_type> *latest, std::unordered_map<coefficient_type, std::pair<bool, bool> *> *isActivePositive,
			std::unordered_map<coefficient_type, List_column *> *columns)
{
    std::list<coefficient_type>::reverse_iterator it;
    std::list<coefficient_type>::iterator it2;
    it = column_->rbegin();
    it++;
    while (it != column_->rend()){
        if (latest->find(*it) != latest->end() && !isActivePositive->at(*it)->first){
	    add(*(columns->at(latest->at(*(it--)))));
        } else if (!isActivePositive->at(*it)->second && !isActivePositive->at(*it)->first) {
            it2 = (++it).base();
            it--; it--;
            column_->erase(it2);
        }
        it++;
    }
}

/**
 * @brief Insert a cell at the end of the column.
 * @param cell value of the cell to be inserted.
 */
inline void List_column::push_back(coefficient_type cell)
{
    column_->push_back(cell);
}

}
}

#endif // LISTCOLUMN_H
