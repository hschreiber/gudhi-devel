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

#ifndef PERSISTENCE_H
#define PERSISTENCE_H

/** @file persistence.h
 * @brief Contains @ref Gudhi::tower_to_filtration::Persistence
 * and @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix classes.
 */

#include <vector>
#include <list>
#include <unordered_map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include <gudhi/tower_converter.h>

namespace Gudhi {
namespace tower_to_filtration {

template<class ComplexStructure, class ColumnType>
/**
 * @brief Takes the elementary operations of a tower in order and computes its persistence barcode.
 */
class Persistence
{
public:
    using vertex = typename ComplexStructure::vertex;	    /**< Type for vertex identifiers. */
    using index = typename Tower_converter<ComplexStructure>::index;	    /**< Type for simplex identifiers. */
    using size_type = typename Tower_converter<ComplexStructure>::size_type; /**< Type for size mesure. */
    using simplex_base = typename std::vector<vertex>;	    /**< Type for simplices. */

    /**
     * @brief Callback function model for output.
     * @param dim dimension of the persistence pair.
     * @param birth birth time of the persistence pair.
     * @param death death time of the persistence pair.
     */
    typedef void (*process_persistence_pair)(int dim, double birth, double death);

    Persistence(size_type reductionInterval, process_persistence_pair outputFunction);
    ~Persistence();

    /**
     * @brief Represents the boundary matrix from which the persistence barcode is computed.
     */
    class Boundary_matrix
    {
    public:
	Boundary_matrix(process_persistence_pair outputFunction);
	~Boundary_matrix();

	void insert_column(index insertionNumber, std::vector<index> &boundary, double timestamp);
	void insert_vertex(index insertionNumber, double timestamp);
	void reduce(size_type start);
	void clear_out();
	void mark_inactive(std::vector<index> &insertionNumbers);
	void mark_inactive(index insertionNumber);

	index get_last_insert_number() const;
	int get_max_dim() const;

    private:
	std::unordered_map<index, ColumnType*> *columns_;			/**< Columns of the matrix. The key is the column number. */
	std::unordered_map<index, index> *latest_;				/**< Pivot to column map. */
	std::unordered_map<index, std::pair<bool, bool>*> *isActivePositive_;	/**< Indicates if a column is active (first value) and/or is positive. */
	std::unordered_map<double, double> *timestamps_;                        /**< Column number to filtration value map. */
	index lastInsertNumber_;						/**< Identifier of the latest inserted simplex (as column). */
	int maxDim_;                                                            /**< Maximal dimension of an inserted simplex. */
	process_persistence_pair outputFunction_;                               /**< Output function. */

	void clear_column(index columnIndex);
	void stream_persistence_pair(int dim, double birth, double death);
    };

    bool add_insertion(simplex_base &simplex, double timestamp);
    bool add_faces_insertions(simplex_base &simplex, double timestamp);
    bool add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim = -1);
    void add_contraction(vertex v, vertex u, double timestamp);
    void finalize_reduction();

private:
    Tower_converter<ComplexStructure> *converter_;  /**< Pointer to @ref Gudhi::tower_to_filtration::Tower_converter<ComplexStructure> */
    Boundary_matrix *matrix_;                       /**< Boundary matrix. */
    size_type reductionInterval_;                   /**< Number of steps between each matrix processing. */
    size_type lastReduction_;                       /**< Number of steps since the last matrix processing. */

    void compute_partial_persistence();
};

template<class ComplexStructure, class ColumnType>
/**
 * @brief Constructor
 * @param reductionInterval number of steps between each matrix processing. The higher the interval, the faster the computing, but the bigger the space consumption.
 * @param outputFunction pointer to a callback function which process the persistence pairs' output stream.
 */
Persistence<ComplexStructure, ColumnType>::Persistence(size_type reductionInterval, process_persistence_pair outputFunction) : reductionInterval_(reductionInterval), lastReduction_(-1)
{
    converter_ = new Tower_converter<ComplexStructure>();
    matrix_ = new Boundary_matrix(outputFunction);
    if (reductionInterval_ < 1) reductionInterval_ = 1;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Destructor
 */
Persistence<ComplexStructure, ColumnType>::~Persistence()
{
    finalize_reduction();
    delete converter_;
    delete matrix_;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Adds an elementary insertion as the next tower operation.
 * @param simplex simplex to be inserted, represented as a vector of its vertex identifiers in increasing order.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @return true if the simplex was not already inserted in the complex, false otherwise.
 */
bool Persistence<ComplexStructure, ColumnType>::add_insertion(simplex_base &simplex, double timestamp)
{
    std::vector<index> boundary;
    index insertionNum;

    if (!converter_->add_insertion(simplex, timestamp, &boundary, &insertionNum)) return false;
    if (simplex.size() == 1) {
        matrix_->insert_vertex(insertionNum, timestamp);
        return true;
    }
    matrix_->insert_column(insertionNum, boundary, timestamp);

    if (fmod(insertionNum, reductionInterval_) == 0) {
        compute_partial_persistence();
    }

    return true;
}

template<class ComplexStructure, class ColumnType>
bool Persistence<ComplexStructure, ColumnType>::add_faces_insertions(simplex_base &simplex, double timestamp)
{
    std::vector<std::vector<index>*> boundaries;
    std::vector<index> insertionNumbers;
    int c = 0;

    bool res = converter_->add_faces_insertions(simplex, timestamp, &boundaries, &insertionNumbers);

    for (std::vector<index>* added : boundaries){
	if (added->size() == 0) {
	    matrix_->insert_vertex(insertionNumbers.at(c), timestamp);
	} else {
	    matrix_->insert_column(insertionNumbers.at(c), *added, timestamp);

	    if (fmod(insertionNumbers.at(c), reductionInterval_) == 0) {
		compute_partial_persistence();
	    }
	}
	++c;
    }

    return res;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Adds a sequence of elementary insertions as the next tower operations. These consists of inserting the edge @p uv (and its vertices if not inserted) and all its possible cofaces.
 * @param u vertex identifier of the first vertex of the edge.
 * @param v vertex identifier of the second vertex of the edge.
 * @param timestamp filtration value for the insertions
 * @param maxExpDim maximal dimension of the cofaces to be inserted ; if -1, then there is no limit.
 * @return true if the edge was not already inserted in the complex, false otherwise.
 */
bool Persistence<ComplexStructure, ColumnType>::add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim)
{
    std::vector<std::vector<index>*> boundaries;
    std::vector<index> insertionNumbers;
    int c = 0;

    bool res = converter_->add_insertions_via_edge_expansion(u, v, timestamp, maxExpDim, &boundaries, &insertionNumbers);

    for (std::vector<index>* boundary : boundaries){
	if (boundary->size() == 0) {
	    matrix_->insert_vertex(insertionNumbers.at(c), timestamp);
	} else {
	    matrix_->insert_column(insertionNumbers.at(c), *boundary, timestamp);

	    if (fmod(insertionNumbers.at(c), reductionInterval_) == 0) {
		compute_partial_persistence();
	    }
	}
	++c;
    }

    return res;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Adds an elementary contraction as the next tower operation.
 * @param v identifier of the contracted vertex which disappears from the complex.
 * @param u identifier of the contracted vertex which remains in the complex.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @exception std::out_of_range If @p v or @p u is not an existing identifier in the current complex ;
 *	Therefore be careful with the order of @p v and @p u to keep coherence with futur contractions.
 */
void Persistence<ComplexStructure, ColumnType>::add_contraction(vertex v, vertex u, double timestamp)
{
    std::vector<std::vector<index>*> boundaries;
    std::vector<index> insertionNumbers;
    bool reduce = false;

    index first = converter_->add_contraction(v, u, timestamp, &boundaries, &insertionNumbers);

    for (typename std::vector<std::vector<index>*>::size_type i = 0; i < boundaries.size(); i++){
	matrix_->insert_column(first + i, *(boundaries.at(i)), timestamp);
        delete boundaries.at(i);
        if (fmod((first + i), reductionInterval_) == 0) reduce = true;
    }
    matrix_->mark_inactive(insertionNumbers);

    if (reduce) compute_partial_persistence();
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Starts matrix processing. To call when the entire tower was given in the input and the last processing has to be made.
 */
inline void Persistence<ComplexStructure, ColumnType>::finalize_reduction()
{
    if (lastReduction_ != matrix_->get_last_insert_number()) matrix_->reduce(lastReduction_ + 1);
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Starts matrix reduction and clearing.
 */
inline void Persistence<ComplexStructure, ColumnType>::compute_partial_persistence()
{
    matrix_->reduce(lastReduction_ + 1);
    matrix_->clear_out();
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Constructor
 * @param outputFunction pointer to a callback function which process the persistence pairs' output stream.
 */
Persistence<ComplexStructure, ColumnType>::Boundary_matrix::Boundary_matrix(process_persistence_pair outputFunction) : lastInsertNumber_(-1), maxDim_(-1), outputFunction_(outputFunction)
{
    columns_ = new std::unordered_map<index, ColumnType*>();
    latest_ = new std::unordered_map<index, index>();
    isActivePositive_ = new std::unordered_map<index, std::pair<bool, bool>*>();
    timestamps_ = new std::unordered_map<double, double>();
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Destructor.
 */
Persistence<ComplexStructure, ColumnType>::Boundary_matrix::~Boundary_matrix()
{
    for (auto it = columns_->begin(); it != columns_->end(); it++){
        delete it->second;
    }
    delete columns_;
    delete latest_;
    for (auto it = isActivePositive_->begin(); it != isActivePositive_->end(); it++){
        delete it->second;
    }
    delete isActivePositive_;
    delete timestamps_;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Insert a nonempty column in the matrix.
 * @param insertionNumber identifier of the simplex whose boundary is represented in this column.
 * @param boundary column to be inserted ; has to represent a new nonzero boundary.
 * @param timestamp filtration value of the corresponding simplex.
 */
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::insert_column(index insertionNumber, std::vector<index> &boundary, double timestamp)
{
    ColumnType *col = new ColumnType(boundary.size() - 1);
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));

    for (int i = 0; i < (int)boundary.size(); i++){
	col->push_back(boundary.at(i));
    }

    columns_->emplace(insertionNumber, col);

    lastInsertNumber_ = insertionNumber;
    if (maxDim_ < (int)boundary.size() - 1) maxDim_ = boundary.size() - 1;
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Insert a vertex in the matrix.
 * @param insertionNumber identifier of the vertex.
 * @param timestamp filtration value of the vertex.
 */
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::insert_vertex(index insertionNumber, double timestamp)
{
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Reduces the matrix from column @p start to the last inserted column in order.
 * @param start number of the column to start the reduction from.
 */
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::reduce(size_type start)
{
    for (int d = maxDim_; d > 0; d--){
	for (size_type i = start; i <= lastInsertNumber_; i++){
            if (columns_->find(i) != columns_->end() && columns_->at(i)->get_dim() == d){
                ColumnType *curr = columns_->at(i);
		index pivot = curr->get_pivot();

                while (pivot != -1 && latest_->find(pivot) != latest_->end()){
		    curr->add(*(columns_->at(latest_->at(pivot))));
                    pivot = curr->get_pivot();
                }

                if (pivot != -1){
                    isActivePositive_->at(i)->second = false;
                    latest_->emplace(pivot, i);
                    clear_column(pivot);
		    stream_persistence_pair(d - 1, pivot, i);
                } else {
                    clear_column(i);
                }
            }
        }
    }
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Clears the matrix from useless cells. See \cite KerberS17 for more details.
 */
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::clear_out()
{
    index r;
    index c;

    for (auto it = columns_->begin(); it != columns_->end(); it++){
        c = it->first;
        ColumnType *column = columns_->at(c);
        r = column->get_pivot();
        if (isActivePositive_->at(r)->first){
            column->clean(latest_, isActivePositive_, columns_);
        }
    }

    for (auto it = columns_->begin(), next_it = columns_->begin(); it != columns_->end(); it = next_it)
    {
        next_it = it; ++next_it;
        c = it->first;
        r = columns_->at(c)->get_pivot();
        if (!isActivePositive_->at(r)->first)
        {
            latest_->erase(r);
            delete isActivePositive_->at(r);
            isActivePositive_->erase(r);
            clear_column(c);
            if (!isActivePositive_->at(c)->first) {
                delete isActivePositive_->at(c);
                isActivePositive_->erase(c);
            }
        }
    }

    for (auto it = isActivePositive_->begin(), next_it = isActivePositive_->begin(); it != isActivePositive_->end(); it = next_it)
    {
        next_it = it; ++next_it;
        c = it->first;
        if (columns_->find(c) == columns_->end() && !it->second->second && !it->second->first)
        {
            delete isActivePositive_->at(c);
            isActivePositive_->erase(c);
        }
    }
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Marks columns as inactive.
 * @param insertionNumbers numbers of the columns to be marked.
 */
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::mark_inactive(std::vector<index> &insertionNumbers)
{
    for (typename std::vector<index>::size_type i = 0; i < insertionNumbers.size(); i++){
	isActivePositive_->at(insertionNumbers.at(i))->first = false;
    }
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Marks a single column as inactive.
 * @param insertionNumber number of the column to be marked.
 */
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::mark_inactive(index insertionNumber)
{
    isActivePositive_->at(insertionNumber)->first = false;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Return the number of the last inserted simplex.
 * @return The number of the last inserted simplex.
 */
inline typename Persistence<ComplexStructure, ColumnType>::index Persistence<ComplexStructure, ColumnType>::Boundary_matrix::get_last_insert_number() const
{
    return lastInsertNumber_;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Returns the maximal dimension of an inserted simplex.
 * @return The maximal dimension of an inserted simplex.
 */
inline int Persistence<ComplexStructure, ColumnType>::Boundary_matrix::get_max_dim() const
{
    return maxDim_;
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Erases the content of a column and the column itself.
 * @param columnIndex number of the column to be deleted.
 */
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::clear_column(index columnIndex)
{
    if (columns_->find(columnIndex) == columns_->end()) return;
    delete columns_->at(columnIndex);
    columns_->erase(columnIndex);
}

template<class ComplexStructure, class ColumnType>
/**
 * @brief Print the given persistence pair (birth and death) into the output.
 * @param dim dimension of the corresponding cycle class.
 * @param birth birth value of the corresponding cycle class.
 * @param death death value of the corresponding cycle class.
 */
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::stream_persistence_pair(int dim, double birth, double death)
{
    if (outputFunction_ != nullptr && timestamps_->at(birth) != timestamps_->at(death)) outputFunction_(dim, timestamps_->at(birth), timestamps_->at(death));
    timestamps_->erase(birth);
    timestamps_->erase(death);
}

}
}

#endif // PERSISTENCE_H
