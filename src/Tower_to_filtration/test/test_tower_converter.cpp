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

#include <iostream>

#define BOOST_TEST_MODULE tower_converter
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/tower_converter.h>
#include <gudhi/hash_complex.h>
#include <gudhi/simplex_tree.h>

using namespace Gudhi::tower_to_filtration;

typedef boost::mpl::list<Hash_complex,Simplex_tree> complex_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(constructor_test, ComplexType, complex_types)
{
    Tower_converter<ComplexType> no_output_tc;

    ComplexType *complex = no_output_tc.get_complex();
    BOOST_CHECK(complex->get_size() == 0);

    BOOST_CHECK(no_output_tc.get_filtration_size() == 0);
    BOOST_CHECK(no_output_tc.get_tower_width() == 0);

    Tower_converter<ComplexType> file_tc("test_tower.txt");

    complex = file_tc.get_complex();
    BOOST_CHECK(complex->get_size() == 0);	//if file not found/open, complex will not be defined and test fails

    BOOST_CHECK(file_tc.get_filtration_size() == 0);
    BOOST_CHECK(file_tc.get_tower_width() == 0);

    std::stringstream ss;
    Tower_converter<ComplexType> stream_tc(&ss);

    complex = stream_tc.get_complex();
    BOOST_CHECK(complex->get_size() == 0);

    BOOST_CHECK(stream_tc.get_filtration_size() == 0);
    BOOST_CHECK(stream_tc.get_tower_width() == 0);
}

template<typename vertex>
bool test_output_stream_first_line(std::stringstream *ss, int *dim, int *filtrationValue, std::vector<vertex> *vertices){
    std::string line;
    if (getline(*ss, line, '\n')){
	std::stringstream nss(line);
	vertices->clear();
	int buf;
	nss >> *dim;
	while (nss >> buf) vertices->push_back(buf);
	BOOST_REQUIRE(!vertices->empty());
	*filtrationValue = vertices->back();
	vertices->pop_back();
	return true;
    } else {
	return false;
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(insertion_test, ComplexType, complex_types)
{
    std::stringstream ss;
    Tower_converter<ComplexType> tc(&ss);
    ComplexType *complex = tc.get_complex();

    std::vector<typename ComplexType::vertex> simplex;
    std::vector<typename ComplexType::index> simplexBoundary;
    typename ComplexType::index simplexInsertionNumber;
    int dim;
    int filtrationValue;
    std::vector<typename ComplexType::vertex> vertices;

    simplex.push_back(0);
    BOOST_CHECK(tc.add_insertion(simplex, 0, &simplexBoundary, &simplexInsertionNumber));
    BOOST_CHECK(simplexBoundary.size() == 0);
    BOOST_CHECK(complex->get_size() == 1);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 0);
    BOOST_CHECK(vertices.size() == 1);
    BOOST_CHECK(filtrationValue == 0);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();     //to enable rewriting in ss.
    BOOST_CHECK(tc.get_filtration_size() == 1);
    BOOST_CHECK(tc.get_tower_width() == 1);

    simplex.at(0) = 1;
    simplexBoundary.clear();
    BOOST_CHECK(tc.add_insertion(simplex, 1, &simplexBoundary, &simplexInsertionNumber));
    BOOST_CHECK(simplexBoundary.size() == 0);
    BOOST_CHECK(complex->get_size() == 2);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 0);
    BOOST_CHECK(vertices.size() == 1);
    BOOST_CHECK(filtrationValue == 1);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();
    BOOST_CHECK(tc.get_filtration_size() == 2);
    BOOST_CHECK(tc.get_tower_width() == 2);

    simplex.at(0) = 0;
    simplex.push_back(1);
    simplexBoundary.clear();
    BOOST_CHECK(tc.add_insertion(simplex, 2, &simplexBoundary, &simplexInsertionNumber));
    BOOST_CHECK(simplexBoundary.size() == 2);
    BOOST_CHECK(simplexBoundary.at(0) == 0);
    BOOST_CHECK(simplexBoundary.at(1) == 1);
    BOOST_CHECK(complex->get_size() == 3);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 1);
    BOOST_CHECK(vertices.size() == 2);
    BOOST_CHECK(filtrationValue == 2);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();
    BOOST_CHECK(tc.get_filtration_size() == 3);
    BOOST_CHECK(tc.get_tower_width() == 3);

    //Simplex already inserted
    simplex.pop_back();
    simplexBoundary.clear();
    BOOST_CHECK(!tc.add_insertion(simplex, 3));

    //Faces output
    ss.str(std::string());
    ss.clear();
    simplexBoundary.clear();
    simplex.clear();
    Tower_converter<ComplexType> tc_faces(&ss, Tower_converter<ComplexType>::FACES);
    complex = tc_faces.get_complex();

    simplex.push_back(0);
    BOOST_CHECK(tc_faces.add_insertion(simplex, 0));
    simplex.at(0) = 1;
    BOOST_CHECK(tc_faces.add_insertion(simplex, 1));
    ss.str(std::string());
    ss.clear();

    simplex.at(0) = 0;
    simplex.push_back(1);
    BOOST_CHECK(tc_faces.add_insertion(simplex, 2, &simplexBoundary, &simplexInsertionNumber));
    BOOST_CHECK(simplexBoundary.size() == 2);
    BOOST_CHECK(simplexBoundary.at(0) == 0);
    BOOST_CHECK(simplexBoundary.at(1) == 1);
    BOOST_CHECK(complex->get_size() == 3);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 1);
    BOOST_CHECK(vertices.size() == 2);
    BOOST_CHECK(filtrationValue == 2);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();
    BOOST_CHECK(tc_faces.get_filtration_size() == 3);
    BOOST_CHECK(tc_faces.get_tower_width() == 3);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(contraction_test, ComplexType, complex_types)
{
    std::stringstream ss;
    Tower_converter<ComplexType> tc(&ss);
    ComplexType *complex = tc.get_complex();

    std::vector<typename ComplexType::vertex> simplex;
    int dim;
    int filtrationValue;
    std::vector<typename ComplexType::vertex> vertices;
    std::vector<std::vector<typename ComplexType::index>*> addedBoundaries;
    std::vector<typename ComplexType::index> removedIndices;

    simplex.push_back(0);
    tc.add_insertion(simplex, 0);
    simplex.at(0) = 1;
    tc.add_insertion(simplex, 1);
    simplex.at(0) = 2;
    tc.add_insertion(simplex, 2);
    simplex.at(0) = 0;
    simplex.push_back(1);
    BOOST_CHECK(tc.add_insertion(simplex, 3));
    simplex.at(0) = 0;
    simplex.at(1) = 2;
    BOOST_CHECK(tc.add_insertion(simplex, 4));
    ss.str(std::string());
    ss.clear();

    typename ComplexType::index first = tc.add_contraction(1, 2, 5, &addedBoundaries, &removedIndices);
    BOOST_CHECK(first != -1);
    BOOST_CHECK(addedBoundaries.size() == 2);
    BOOST_CHECK(removedIndices.size() == 4);
    BOOST_CHECK(complex->get_size() == 3);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 1);
    BOOST_CHECK(vertices.size() == 2);
    BOOST_CHECK(filtrationValue == 5);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 2);
    BOOST_CHECK(vertices.size() == 3);
    BOOST_CHECK(filtrationValue == 5);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();
    BOOST_CHECK(tc.get_filtration_size() == 7);
    BOOST_CHECK(tc.get_tower_width() == 5);

    //Faces output
    ss.str(std::string());
    ss.clear();
    simplex.clear();
    for (std::vector<typename ComplexType::index>* v : addedBoundaries) delete v;
    addedBoundaries.clear();
    removedIndices.clear();
    Tower_converter<ComplexType> tc_faces(&ss, Tower_converter<ComplexType>::FACES);
    complex = tc_faces.get_complex();

    simplex.push_back(0);
    tc_faces.add_insertion(simplex, 0);
    simplex.at(0) = 1;
    tc_faces.add_insertion(simplex, 1);
    simplex.at(0) = 2;
    tc_faces.add_insertion(simplex, 2);
    simplex.at(0) = 0;
    simplex.push_back(1);
    BOOST_CHECK(tc_faces.add_insertion(simplex, 3));
    simplex.at(0) = 0;
    simplex.at(1) = 2;
    BOOST_CHECK(tc_faces.add_insertion(simplex, 4));
    ss.str(std::string());
    ss.clear();

    first = tc_faces.add_contraction(1, 2, 5, &addedBoundaries, &removedIndices);
    BOOST_CHECK(first != -1);
    BOOST_CHECK(addedBoundaries.size() == 2);
    BOOST_CHECK(removedIndices.size() == 4);
    BOOST_CHECK(complex->get_size() == 3);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 1);
    BOOST_CHECK(vertices.size() == 2);
    BOOST_CHECK(filtrationValue == 5);
    test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices);
    BOOST_CHECK(dim == 2);
    BOOST_CHECK(vertices.size() == 3);
    BOOST_CHECK(filtrationValue == 5);
    BOOST_CHECK(!test_output_stream_first_line(&ss, &dim, &filtrationValue, &vertices));
    ss.clear();
    BOOST_CHECK(tc_faces.get_filtration_size() == 7);
    BOOST_CHECK(tc_faces.get_tower_width() == 5);
}

