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

#ifndef CONCEPT_COMPLEX_STRUCTURE_H_
#define CONCEPT_COMPLEX_STRUCTURE_H_

/** @file ComplexStructure.h
 * @brief Contains @ref Gudhi::tower_to_filtration::ComplexStructure concept.
 */

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Data structure for a complex. This concept describes the methodes required by the module.
 */
class ComplexStructure{
public:
    typename vertex;		    /**< Type for vertex identifiers. Should be able to be used by std::unordered_map as key type without customized hash function. Needs to be comparable and have a total order (<, ==). */
    typename simplex_handle;	    /**< Type for simplex identifiers. Should be able to be used by std::unordered_map as key type without customized hash function. */
    typename simplex_vertex_range;  /**< Ordered range of vertex identifiers. Should be iteratable with a begin() and end() function. */

    /**
     * @brief Constructor with no parameters
     */
    ComplexStructure();

    /**
     * @brief Insert a simplex in the complex, whose facets were already inserted in the complex.
     * @param simplex simplex to be inserted, described as the vector of its vertex identifiers. The module will give them in increasing order.
     * @param handle pointer to a simplex identifier ; the method should store here the identifier of the just inserted simplex.
     *	    This identifier needs to be unique and fixed: it should not change later in the process and a same simplex should always be associated with the same identifier.
     *	    If the simplex is removed from the simplex, the identifier can be invalidated.
     * @return true if the simplex was not already inserted in the complex, false otherwise.
     *	    This module uses the result in an `if` statement, the output does not have to be `bool` if the test result is the same.
     */
    bool insert_simplex(std::vector<vertex> &simplex, simplex_handle *handle);
    /**
     * @brief Insert a simplex in the complex and all its facets if they are not already inserted.
     * @param simplex simplex simplex to be inserted, described as the vector of its vertex identifiers. The module will give them in increasing order.
     * @param addedSimplices pointer to an (empty) vector of simplex identifiers ; the method should store in a valid filtration order all newly inserted simplices here.
     *	    These identifiers need to be unique and fixed: they should not change later in the process and a same simplex should always be associated with the same identifier.
     *	    If the simplex is removed from the simplex, the identifier can be invalidated.
     * @return true if @p simplex was not already inserted in the complex, false otherwise.
     */
    bool insert_simplex_and_faces(std::vector<vertex> &simplex, std::vector<simplex_handle> *addedSimplices);
    /**
     * @brief Inserts the edge @p uv (and its vertices if not inserted) and all its possible cofaces.
     * @param u vertex identifier of the first vertex of the edge ; the module calls the method such that @p u < @p v.
     * @param v vertex identifier of the second vertex of the edge ; the module calls the method such that @p u < @p v.
     * @param maxDim the maximal dimension of the cofaces to be inserted ; if -1, then there is no limit.
     * @param addedSimplices pointer to an (empty) vector of simplex identifiers ; the method should store in a valid filtration order all newly inserted simplices here.
     *	    These identifiers need to be unique and fixed: they should not change later in the process and a same simplex should always be associated with the same identifier.
     *	    If the simplex is removed from the simplex, the identifier can be invalidated.
     * @return true if the edge was not already inserted in the complex, false otherwise.
     */
    bool insert_edge_and_expand(vertex u, vertex v, int maxDim, std::vector<simplex_handle> *addedSimplices);
    /**
     * @brief Remove a simplex and all its cofaces from the complex.
     * @param vertex simplex to be removed. This modules only needs to call the function on vertices.
     * @param removedIndices pointer to an (empty) vector of simplex identifiers ; the identifiers of all the removed simplices should be stored in the vector.
     */
    void remove_simplex(simplex_handle &vertex, std::vector<simplex_handle> *removedIndices);
    /**
     * @brief Compute the smallest closed star of the vertices v and u and returns the corresponding vertex
     * in a time depending on the smallest closed star and not the biggest one.
     *
     * The star of a simplex is the set of simplices containing this simplex. The closed star is the smallest subcomplex containing the star.
     *
     * @param v identifier of the first vertex
     * @param u identifier of the second vertex
     * @param closedStar pointer to an empty vector of simplices ; the method should store here the simplices contained in the smallest closed star,
     *	    ordered such that possible faces of a simplex come before the simplex itself.
     * @return The identifier of the vertex which had the smallest closed star.
     */
    simplex_handle get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> *closedStar);
    /**
     * @brief Computes the boundary of a simplex and returns the identifier of the simplex.
     * @param simplex identifier of the simplex which boundary will be computed.
     * @param boundary pointer to an (empty) vector of simplex identifiers ; the method should store here the simplices in the boundary in increasing order of insertion.
     * @return The identifier of the simplex.
     */
    void get_boundary(simplex_handle &simplex, std::vector<simplex_handle> *boundary);
    /**
     * @brief Returns a range of the identifiers of the given simplex' vertices.
     * @param simplex identifier of the simplex of which the vertices should be returned.
     * @return A range of the identifiers of the given simplex' vertices.
     */
    simplex_vertex_range get_vertices(simplex_handle &simplex);
};

}
}

#endif  // CONCEPT_COMPLEX_STRUCTURE_H_
