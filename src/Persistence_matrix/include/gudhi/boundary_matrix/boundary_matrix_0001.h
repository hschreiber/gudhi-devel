/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BOUNDARY_MATRIX_0001_H
#define BOUNDARY_MATRIX_0001_H

#include <iostream>
#include <vector>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Boundary_matrix_with_removals
		: public Master_matrix::Base_swap_option,
		  public Master_matrix::Base_pairing_option
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using Row_type = void;

	Boundary_matrix_with_removals();
	template<class Boundary_type = boundary_type>
	Boundary_matrix_with_removals(const std::vector<Boundary_type>& orderedBoundaries);
	Boundary_matrix_with_removals(unsigned int numberOfColumns);
	Boundary_matrix_with_removals(const Boundary_matrix_with_removals& matrixToCopy);
	Boundary_matrix_with_removals(Boundary_matrix_with_removals&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type get_row(index rowIndex) const;
	void remove_maximal_simplex(index columnIndex);		//does not update barcode!

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex);

	Boundary_matrix_with_removals& operator=(Boundary_matrix_with_removals other);
	friend void swap(Boundary_matrix_with_removals& matrix1,
					 Boundary_matrix_with_removals& matrix2){
		swap(static_cast<typename Master_matrix::Base_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Base_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_pairing_option&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		matrix1.dimensions_.swap(matrix2.dimensions_);
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Base_swap_option;
	using pair_opt = typename Master_matrix::Base_pairing_option;

	typename Master_matrix::column_container_type matrix_;
	std::vector<unsigned int> dimensions_;
	dimension_type maxDim_;
	index nextInsertIndex_;
};

template<class Master_matrix>
inline Boundary_matrix_with_removals<Master_matrix>::Boundary_matrix_with_removals()
	: Master_matrix::Base_swap_option(matrix_),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Boundary_matrix_with_removals<Master_matrix>::Boundary_matrix_with_removals(const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Base_swap_option(matrix_, orderedBoundaries.size()),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  matrix_(orderedBoundaries.size()),
	  nextInsertIndex_(orderedBoundaries.size())
{
	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		const Boundary_type& b = orderedBoundaries[i];
		matrix_.try_emplace(i, Column_type(b));

		unsigned int dim = (b.size() == 0) ? 0 : static_cast<int>(b.size()) - 1;
		if (dimensions_.size() <= dim) dimensions_.resize(dim + 1);
		++(dimensions_[dim]);
	}

	maxDim_ = dimensions_.size() - 1;
}

template<class Master_matrix>
inline Boundary_matrix_with_removals<Master_matrix>::Boundary_matrix_with_removals(unsigned int numberOfColumns)
	: Master_matrix::Base_swap_option(matrix_, numberOfColumns),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  matrix_(numberOfColumns),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
inline Boundary_matrix_with_removals<Master_matrix>::Boundary_matrix_with_removals(const Boundary_matrix_with_removals &matrixToCopy)
	: Master_matrix::Base_swap_option(matrixToCopy),
	  Master_matrix::Base_pairing_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  dimensions_(matrixToCopy.dimensions_),
	  maxDim_(matrixToCopy.maxDim_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
	if constexpr (pair_opt::isActive_){
		pair_opt::matrix_ = &matrix_;
		pair_opt::maxDim_ = &maxDim_;
	}
}

template<class Master_matrix>
inline Boundary_matrix_with_removals<Master_matrix>::Boundary_matrix_with_removals(Boundary_matrix_with_removals &&other) noexcept
	: Master_matrix::Base_swap_option(std::move(other)),
	  Master_matrix::Base_pairing_option(std::move(other)),
	  matrix_(std::move(other.matrix_)),
	  dimensions_(std::move(other.dimensions_)),
	  maxDim_(std::exchange(other.maxDim_,-1)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
	if constexpr (pair_opt::isActive_){
		pair_opt::matrix_ = &matrix_;
		pair_opt::maxDim_ = &maxDim_;
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void Boundary_matrix_with_removals<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	if constexpr (swap_opt::isActive_){
		swap_opt::indexToRow_[nextInsertIndex_] = nextInsertIndex_;
		swap_opt::rowToIndex_[nextInsertIndex_] = nextInsertIndex_;
	}

	matrix_.try_emplace(nextInsertIndex_++, boundary);

	unsigned int dim = (boundary.size() == 0) ? 0 : boundary.size() - 1;
	if (dimensions_.size() <= dim) dimensions_.resize(dim + 1);
	++(dimensions_[dim]);
	maxDim_ = dimensions_.size() - 1;
}

template<class Master_matrix>
inline typename Boundary_matrix_with_removals<Master_matrix>::Column_type &Boundary_matrix_with_removals<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline const typename Boundary_matrix_with_removals<Master_matrix>::Column_type &Boundary_matrix_with_removals<Master_matrix>::get_column(index columnIndex) const
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Boundary_matrix_with_removals<Master_matrix>::Row_type Boundary_matrix_with_removals<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::remove_maximal_simplex(index columnIndex)
{
	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

	int dim = matrix_.at(columnIndex).get_dimension();
	--(dimensions_[dim]);
	while (dimensions_.back() == 0)
		dimensions_.pop_back();
	maxDim_ = dimensions_.size() - 1;

	matrix_.erase(columnIndex);

	if constexpr (swap_opt::isActive_){
		auto it = swap_opt::indexToRow_.find(columnIndex);
		swap_opt::rowToIndex_.erase(it->second);
		swap_opt::indexToRow_.erase(it);
	}
	if constexpr (pair_opt::isActive_){
		if (pair_opt::isReduced_){
			auto bar = pair_opt::indexToBar_.at(columnIndex);

			if (bar->death == -1) pair_opt::barcode_.erase(bar);
			else bar->death = -1;

			pair_opt::indexToBar_.erase(columnIndex);
		}
	}
}

template<class Master_matrix>
inline dimension_type Boundary_matrix_with_removals<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Boundary_matrix_with_removals<Master_matrix>::get_number_of_columns() const
{
	return nextInsertIndex_;
}

template<class Master_matrix>
inline dimension_type Boundary_matrix_with_removals<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_.at(columnIndex).get_dimension();
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += sourceColumn;
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		matrix_.at(columnIndex).clear(swap_opt::indexToRow_.at(rowIndex));
	} else {
		matrix_.at(columnIndex).clear(rowIndex);
	}
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::zero_column(index columnIndex)
{
	matrix_[columnIndex].clear();
}

template<class Master_matrix>
inline bool Boundary_matrix_with_removals<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		return !(matrix_.at(columnIndex).is_non_zero(swap_opt::indexToRow_.at(rowIndex)));
	} else {
		return !(matrix_.at(columnIndex).is_non_zero(rowIndex));
	}
}

template<class Master_matrix>
inline bool Boundary_matrix_with_removals<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_.at(columnIndex).is_empty();
}

template<class Master_matrix>
inline index Boundary_matrix_with_removals<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_column_with_pivot' is not implemented for the chosen options.");
	return 0;
}

template<class Master_matrix>
inline index Boundary_matrix_with_removals<Master_matrix>::get_pivot(index columnIndex)
{
	return matrix_.at(columnIndex).get_pivot();
}

template<class Master_matrix>
inline Boundary_matrix_with_removals<Master_matrix> &Boundary_matrix_with_removals<Master_matrix>::operator=(Boundary_matrix_with_removals other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	matrix_.swap(other.matrix_);
	dimensions_.swap(other.dimensions_);
	std::swap(maxDim_, other.maxDim_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void Boundary_matrix_with_removals<Master_matrix>::print()
{
	std::cout << "Boundary_matrix_with_removals:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = matrix_.at(i);
		for (const auto e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BOUNDARY_MATRIX_0001_H