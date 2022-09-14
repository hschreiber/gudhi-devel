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

#ifndef TOWER_CONVERTER_H
#define TOWER_CONVERTER_H

/** @file tower_converter.h
 * @brief Contains @ref Gudhi::tower_to_filtration::Tower_converter class.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

namespace Gudhi {
namespace tower_to_filtration {

template<class ComplexStructure>
/**
 * @brief Takes the elementary operations of a tower in order and convert them into an equivalent filtration.
 */
class Tower_converter
{
public:
	using vertex = typename ComplexStructure::vertex;			/**< Type for vertex identifiers. */
	using simplex_handle = typename ComplexStructure::simplex_handle;
	using simplex_vertex_range = typename ComplexStructure::simplex_vertex_range;
	using index = long long;									/**< Type for simplex identifiers. */
	using size_type = long long;								/**< Type for size mesure. */
	using simplex_base = typename std::vector<vertex>;			/**< Type for simplices. */

	/**
	 * @brief Callback function model for output.
	 * @param complex pointer to the stored complex.
	 * @param simplex vector of vertex identifiers of the next simplex in the filtration.
	 * @param timestamp filtration value of the next simplex in the filtration.
	 */
	typedef void (*process_output)(ComplexStructure &complex, simplex_handle &simplex, double timestamp);

	Tower_converter();
	Tower_converter(process_output outputFunction);

	bool add_insertion(simplex_base &simplex, double timestamp);
	bool add_insertion(simplex_base &simplex, double timestamp,
					   std::vector<index> &simplexBoundary, index &simplexInsertionNumber);
	bool add_faces_insertions(simplex_base &simplex, double timestamp);
	bool add_faces_insertions(simplex_base &simplex, double timestamp,
							  std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers);
	bool add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim = -1);
	bool add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim,
										   std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers);
	index add_contraction(vertex v, vertex u, double timestamp);
	index add_contraction(vertex v, vertex u, double timestamp,
						  std::vector<std::vector<index> > &addedBoundaries, std::vector<index> &removedIndices);

	const ComplexStructure& get_complex() const;
	size_type get_filtration_size() const;
	size_type get_tower_width() const;

private:
	ComplexStructure complex_;						/**< Current complex. */
	std::unordered_map<vertex, vertex> vertices_;	/**< Current vertices in the complex. Keeps the coherence between vertex identifiers outside and inside the class. */
	process_output outputFunction_;					/**< Output function. */
	size_type filtrationSize_;						/**< Current filtration size. */
	size_type towerWidth_;							/**< Current tower width. */
	size_type complexSize_;
	index maxIndex_;

	simplex_base _translate_simplex(simplex_base &simplex);
	void _add_cone(std::vector<simplex_handle> &simplices, double timestamp);
	void _add_cone(std::vector<simplex_handle> &simplices, double timestamp, std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers);
	void _add_cone(simplex_handle &simplex, double timestamp, std::vector<index> &boundary, index &insertionNumber);
	void _get_union(vertex v, std::vector<simplex_handle> &simplices, std::vector<simplex_base> &unions);
	void _stream_simplex(simplex_handle &simplex, double timestamp);
};

template<class ComplexStructure>
/**
 * @brief Constructor without parameters.
 *
 * Initializes the members. The output option is set as "no output".
 */
Tower_converter<ComplexStructure>::Tower_converter()
	: outputFunction_(nullptr), filtrationSize_(0), towerWidth_(0), complexSize_(0), maxIndex_(-1)
{}

template<class ComplexStructure>
/**
 * @brief Full constructor.
 *
 * Initializes the members. The output stream will be redirected to @p outputFunction.
 *
 * @param outputFunction pointer to a callback function which process the filtration output stream.
 */
Tower_converter<ComplexStructure>::Tower_converter(process_output outputFunction)
	: outputFunction_(outputFunction), filtrationSize_(0), towerWidth_(0), complexSize_(0), maxIndex_(-1)
{}

template<class ComplexStructure>
/**
 * @brief Adds an elementary insertion as the next tower operation and convert it into the output stream.
 * @param simplex simplex to be inserted, represented as a vector of its vertex identifiers in increasing order.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @return true if the simplex was not already inserted in the complex, false otherwise.
 */
inline bool Tower_converter<ComplexStructure>::add_insertion(simplex_base &simplex, double timestamp)
{
	simplex_base transSimplex = _translate_simplex(simplex);
	std::vector<simplex_handle> handle = complex_.insert_simplex(transSimplex);

	if (handle.empty()) return false;

	_add_cone(handle, timestamp);

	return true;
}

template<class ComplexStructure>
/**
 * @brief Adds an elementary insertion as the next tower operation and convert it into the output stream.
 * @param simplex simplex to be inserted, represented as a vector of its vertex identifiers in increasing order.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @param simplexBoundary pointer to an (empty) vector, where the identifiers of the boundary of the inserted simplex will be stored.
 * @param simplexInsertionNumber pointer to an identifier, which will be replaced by the one of the inserted simplex.
 * @return true if the simplex was not already inserted in the complex, false otherwise.
 */
bool Tower_converter<ComplexStructure>::add_insertion(simplex_base &simplex, double timestamp, std::vector<index> &simplexBoundary, index &simplexInsertionNumber)
{
	simplex_base transSimplex = _translate_simplex(simplex);
	std::vector<simplex_handle> handle = complex_.insert_simplex(transSimplex);

	if (handle.empty()) return false;

	_add_cone(handle[0], timestamp, simplexBoundary, simplexInsertionNumber);
	if (complexSize_ > towerWidth_) towerWidth_ = complexSize_;

	return true;
}

template<class ComplexStructure>
bool Tower_converter<ComplexStructure>::add_faces_insertions(simplex_base &simplex, double timestamp)
{
	simplex_base transSimplex = _translate_simplex(simplex);
	std::vector<simplex_handle> insertedSimplices = complex_.insert_simplex_and_faces(transSimplex);

	_add_cone(insertedSimplices, timestamp);

	return !insertedSimplices.empty();
}

template<class ComplexStructure>
bool Tower_converter<ComplexStructure>::add_faces_insertions(simplex_base &simplex, double timestamp, std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers)
{
	simplex_base transSimplex = _translate_simplex(simplex);
	std::vector<simplex_handle> insertedSimplices = complex_.insert_simplex_and_faces(transSimplex);

	_add_cone(insertedSimplices, timestamp, boundaries, insertionNumbers);

	return !insertedSimplices.empty();
}

template<class ComplexStructure>
/**
 * @brief Adds a sequence of elementary insertions as the next tower operations. These consists of inserting the edge @p uv (and its vertices if not inserted) and all its possible cofaces.
 * @param u vertex identifier of the first vertex of the edge.
 * @param v vertex identifier of the second vertex of the edge.
 * @param timestamp filtration value for the insertions
 * @param maxExpDim maximal dimension of the cofaces to be inserted ; if -1, then there is no limit.
 * @return true if the edge was not already inserted in the complex, false otherwise.
 */
bool Tower_converter<ComplexStructure>::add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim)
{
	vertex first;
	vertex second;

	vertices_.emplace(u, u);
	vertices_.emplace(v, v);

	first = vertices_[u] < vertices_[v] ? vertices_[u] : vertices_[v];
	second = vertices_[u] < vertices_[v] ? vertices_[v] : vertices_[u];

	std::vector<simplex_handle> insertedSimplices = complex_.insert_edge_and_expand(first, second, maxExpDim);

	_add_cone(insertedSimplices, timestamp);

	return !insertedSimplices.empty();
}

template<class ComplexStructure>
/**
 * @brief Adds a sequence of elementary insertions as the next tower operations. These consists of inserting the edge @p uv (and its vertices if not inserted) and all its possible cofaces.
 * @param u vertex identifier of the first vertex of the edge.
 * @param v vertex identifier of the second vertex of the edge.
 * @param timestamp filtration value for the insertions
 * @param maxExpDim maximal dimension of the cofaces to be inserted ; if -1, then there is no limit.
 * @param addedSimplices pointer to an (empty) vector of simplices ; the method stores in insertion order all newly inserted simplices here.
 * @param boundaries pointer to an (empty) vector of boundary identifiers ; the method stores there the identifiers of the
 *	boundary simplices of the simplices in @p addedSimplices in the corresponding order.
 * @param insertionNumbers pointer to an (empty) vector of identifiers ;
 *	the method stores there the identifiers of the simplices in @p addedSimplices in the corresponding order.
 * @return
 */
bool Tower_converter<ComplexStructure>::add_insertions_via_edge_expansion(vertex u, vertex v, double timestamp, int maxExpDim,
																		  std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers)
{
	vertex first;
	vertex second;

	vertices_.emplace(u, u);
	vertices_.emplace(v, v);

	first = vertices_[u] < vertices_[v] ? vertices_[u] : vertices_[v];
	second = vertices_[u] < vertices_[v] ? vertices_[v] : vertices_[u];

	std::vector<simplex_handle> insertedSimplices = complex_.insert_edge_and_expand(first, second, maxExpDim);

	_add_cone(insertedSimplices, timestamp, boundaries, insertionNumbers);

	return !insertedSimplices.empty();
}

template<class ComplexStructure>
/**
 * @brief Adds an elementary contraction as the next tower operation and convert it into a equivalent sequence of insertions into the output stream.
 * @param v identifier of the contracted vertex which disappears from the complex.
 * @param u identifier of the contracted vertex which remains in the complex.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @exception std::out_of_range If @p v or @p u is not an existing identifier in the current complex ;
 *	Therefore be careful with the order of @p v and @p u to keep coherence with futur contractions.
 * @return The identifier of the first new simplex in the equivalent insertion ; the remaining new simplices will take the identifiers which follows continuously.
 */
inline typename Tower_converter<ComplexStructure>::index Tower_converter<ComplexStructure>::add_contraction(vertex v, vertex u, double timestamp)
{
	std::vector<simplex_handle> closedStar;
	std::vector<simplex_base> unions;
	vertex tv = vertices_.at(v), tu = vertices_.at(u);
	simplex_handle disappearing = complex_.get_smallest_closed_star(tv, tu, closedStar);
	index first = -1;

	vertices_.erase(v);
	if (*(complex_.get_vertices(disappearing).begin()) == tu){
		vertices_.at(u) = tv;
		_get_union(tv, closedStar, unions);
	} else {
		_get_union(tu, closedStar, unions);
	}

	for (auto it = unions.begin(); it != unions.end(); it++){
		std::vector<simplex_handle> handle = complex_.insert_simplex(*it);
		if (!handle.empty()) {
			_stream_simplex(handle[0], timestamp);
			if (first == -1) first = maxIndex_;
		}
	}

	complexSize_ -= complex_.remove_simplex(disappearing).size();

	if (complexSize_ > towerWidth_) towerWidth_ = complexSize_;

	return first;
}

template<class ComplexStructure>
/**
 * @brief Adds an elementary contraction as the next tower operation and convert it into a equivalent sequence of insertions into the output stream.
 * @param v identifier of the contracted vertex which disappears from the complex.
 * @param u identifier of the contracted vertex which remains in the complex.
 * @param timestamp time value or filtration value which will be associated to the operation in the filtration. Has to be equal or higher to the precedent ones.
 * @param addedBoundaries pointer to an (empty) vector, where the boundaries of the inserted simplices will be stored.
 * @param removedIndices pointer to an (empty) vector, where the identifiers of the simplices which become inactive will be stored.
 * @exception std::out_of_range If @p v or @p u is not an existing identifier in the current complex ;
 *	Therefore be careful with the order of @p v and @p u to keep coherence with futur contractions.
 * @return The identifier of the first new simplex in the equivalent insertion ; the remaining new simplices will take the identifiers which follows continuously.
 */
typename Tower_converter<ComplexStructure>::index Tower_converter<ComplexStructure>::add_contraction(
		vertex v, vertex u, double timestamp,
		std::vector<std::vector<index> > &addedBoundaries, std::vector<index> &removedIndices)
{
	std::vector<simplex_handle> closedStar;
	std::vector<simplex_base> unions;
	vertex tv = vertices_.at(v), tu = vertices_.at(u);
	simplex_handle disappearing = complex_.get_smallest_closed_star(tv, tu, closedStar);
	index first = -1;

	vertices_.erase(v);
	if (*(complex_.get_vertices(disappearing).begin()) == tu){
		vertices_.at(u) = tv;
		_get_union(tv, closedStar, unions);
	} else {
		_get_union(tu, closedStar, unions);
	}

	for (auto it = unions.begin(); it != unions.end(); it++){
		std::vector<simplex_handle> handle = complex_.insert_simplex(*it);
		if (!handle.empty()) {
			index insertionNumber;
			addedBoundaries.push_back(std::vector<index>());
			_add_cone(handle[0], timestamp, addedBoundaries.back(), insertionNumber);

			if (first == -1) first = insertionNumber;
		}
	}

	removedIndices = complex_.remove_simplex(disappearing);
	complexSize_ -= removedIndices.size();

	if (complexSize_ > towerWidth_) towerWidth_ = complexSize_;

	return first;
}

template<class ComplexStructure>
/**
 * @brief Returns the current size of the filtration.
 * @return The current size of the filtration.
 */
inline typename Tower_converter<ComplexStructure>::size_type Tower_converter<ComplexStructure>::get_filtration_size() const
{
	return filtrationSize_;
}

template<class ComplexStructure>
/**
 * @brief Returns the maximal size of the complex until now.
 * @return The maximal size of the complex until now.
 */
inline typename Tower_converter<ComplexStructure>::size_type Tower_converter<ComplexStructure>::get_tower_width() const
{
	return towerWidth_;
}

template<class ComplexStructure>
inline typename Tower_converter<ComplexStructure>::simplex_base Tower_converter<ComplexStructure>::_translate_simplex(simplex_base &simplex)
{
	simplex_base transSimplex(simplex.size());

	unsigned int i = 0;
	for (vertex v : simplex){
		if (vertices_.find(v) == vertices_.end()){
			transSimplex[i++] = v;
			vertices_.emplace(v, v);
		} else {
			transSimplex[i++] = vertices_.at(v);
		}
	}
	std::sort(transSimplex.begin(), transSimplex.end());

	return transSimplex;
}

template<class ComplexStructure>
inline void Tower_converter<ComplexStructure>::_add_cone(std::vector<simplex_handle> &simplices, double timestamp)
{
	for (simplex_handle added : simplices){
		_stream_simplex(added, timestamp);
	}

	if (complexSize_ > towerWidth_) towerWidth_ = complexSize_;
}

template<class ComplexStructure>
inline void Tower_converter<ComplexStructure>::_add_cone(std::vector<simplex_handle> &simplices, double timestamp, std::vector<std::vector<index> > &boundaries, std::vector<index> &insertionNumbers)
{
	boundaries.resize(simplices.size());
	insertionNumbers.resize(simplices.size());

	for (unsigned int i = 0; i < simplices.size(); ++i){
		_add_cone(simplices[i], timestamp, boundaries[i], insertionNumbers[i]);
	}

	if (complexSize_ > towerWidth_) towerWidth_ = complexSize_;
}

template<class ComplexStructure>
inline void Tower_converter<ComplexStructure>::_add_cone(simplex_handle &simplex, double timestamp, std::vector<index> &boundary, index &insertionNumber)
{
	_stream_simplex(simplex, timestamp);

	insertionNumber = maxIndex_;
	std::vector<simplex_handle> trans_boundary = complex_.get_boundary(simplex);
	boundary.resize(trans_boundary.size());
	unsigned int j = 0;
	for (simplex_handle h : trans_boundary) {
		boundary[j++] = complex_.get_key(h);
	}
}

template<class ComplexStructure>
/**
 * @brief Return pointer to the stored complex.
 * @return The pointer to the stored complex.
 */
inline const ComplexStructure& Tower_converter<ComplexStructure>::get_complex() const
{
	return complex_;
}

template<class ComplexStructure>
/**
 * @brief Make the union of a set of simplices and a vertex.
 * @param v vertex identifier to unify.
 * @param simplices vector of simplices to unify with @p v.
 * @param unions pointer to an empty vector of simplices. The method will fill the vector with the resulting simplices.
 */
void Tower_converter<ComplexStructure>::_get_union(vertex v, std::vector<simplex_handle> &simplices, std::vector<simplex_base> &unions)
{
	unions.resize(simplices.size());
	int c = 0;

	for (simplex_handle handle : simplices){
		simplex_vertex_range vertices = complex_.get_vertices(handle);
		auto it = vertices.begin();
		while (it != vertices.end() && *it < v){
			unions.at(c).push_back(*it);
			++it;
		}
		if ((it != vertices.end() && *it != v) || it == vertices.end()) unions.at(c).push_back(v);
		while (it != vertices.end()){
			unions.at(c).push_back(*it);
			++it;
		}
		++c;
	}
}

template<class ComplexStructure>
/**
 * @brief Writes the simplex as an insertion in the output stream.
 * @param simplex simplex to be inserted.
 * @param timestamp filtration value of the insertion.
 */
void Tower_converter<ComplexStructure>::_stream_simplex(simplex_handle &simplex, double timestamp)
{
	++maxIndex_;
	++complexSize_;
	++filtrationSize_;
	if (outputFunction_ != nullptr) outputFunction_(complex_, simplex, timestamp);
}

}
}

#endif // TOWER_CONVERTER_H
