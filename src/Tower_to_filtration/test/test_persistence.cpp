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
#include <string>

#define BOOST_TEST_MODULE persistence
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Persistent_homology/persistence.h>
#include <gudhi/Persistent_homology/heap_column.h>
#include <gudhi/Persistent_homology/list_column.h>
#include <gudhi/hash_complex.h>
#include <gudhi/simplex_tree.h>
#include <gudhi/tc_reading_utilities.h>

using namespace Gudhi::tower_to_filtration;

typedef boost::mpl::list<Hash_complex,Simplex_tree> complex_types;

// For real unitary tests, access to matrix is needed to verify each call to add_insertion or add_contraction. TODO? (friend class ? public matrix ?)

bool verify_resulting_file(){
    std::ifstream file("res.txt");
    std::string line;
    int dim;
    int lineNumber = 0;

    if (file.is_open()){
	while (getline(file, line, '\n')){
	    std::stringstream ss(line);
	    lineNumber++;
	    ss >> dim;
	    if (dim != 0) return false;
	}
	file.close();
    } else {
	std::cout << "Unable to open input file\n";
	file.setstate(std::ios::failbit);
	return false;
    }

    return true;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(list_test, ComplexType, complex_types)
{
    Persistence<ComplexType, List_column> pers(5, "res.txt");
    std::ifstream file("test_tower.txt");

    if (file.is_open()){
	file >> pers;
	file.close();
    } else {
	std::cout << "Unable to open input file\n";
	file.setstate(std::ios::failbit);
    }

    BOOST_CHECK(verify_resulting_file());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(heap_test, ComplexType, complex_types)
{
    Persistence<ComplexType, Heap_column> pers(5, "res.txt");
    std::ifstream file("test_tower.txt");

    if (file.is_open()){
	file >> pers;
	file.close();
    } else {
	std::cout << "Unable to open input file\n";
	file.setstate(std::ios::failbit);
    }

    BOOST_CHECK(verify_resulting_file());
}

