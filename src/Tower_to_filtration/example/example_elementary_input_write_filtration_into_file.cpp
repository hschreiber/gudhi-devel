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
#include "gudhi/tower_converter.h"
#include "gudhi/hash_complex.h"
#include "gudhi/tc_reading_utilities.h"

using namespace Gudhi::tower_to_filtration;	//module namespace

/**
 * @brief Print usage of example file
 */
void print_usage(){
	std::cout << "Usage:" << std::endl;
	std::cout << "  ./example_elementary_input_write_filtration_into_file output_file_name" << std::endl;
}

using vertex = Hash_complex::vertex;

std::ofstream *outputFile;		/**< Pointer to output file. */

void write_tc_results_to_file_as_vertices(Hash_complex *complex, std::vector<vertex> &simplex, double timestamp);	/**< Example 1 of callback function for Tower_converter. */
void write_tc_results_to_file_as_faces(Hash_complex *complex, std::vector<vertex> &simplex, double timestamp);		/**< Example 2 of callback function for Tower_converter. */

int main(int argc, char *argv[])
{
	if (argc != 2){
		print_usage();
		return 0;
	}

	std::vector<vertex> vertices;
	outputFile = new std::ofstream(argv[1]);

	if (!outputFile->is_open()){
		std::cout << "Unable to open output file." << std::endl;
		return 0;
	}

	Tower_converter<Hash_complex> tc(write_tc_results_to_file_as_vertices);	// output callback function ; another example: write_tc_results_to_file_as_faces.

	vertices.push_back(0);
	tc.add_insertion(vertices, 0);		// add vertex 0 at time 0

	vertices.at(0) = 1;
	tc.add_insertion(vertices, 1);		// add vertex 1 at time 1

	vertices.at(0) = 2;
	tc.add_insertion(vertices, 2);		// add vertex 2 at time 2

	vertices.at(0) = 3;
	tc.add_insertion(vertices, 3);		// add vertex 3 at time 3

	vertices.at(0) = 4;
	tc.add_insertion(vertices, 4);		// add vertex 4 at time 4

	vertices.at(0) = 5;
	tc.add_insertion(vertices, 5);		// add vertex 5 at time 5

	vertices.at(0) = 6;
	tc.add_insertion(vertices, 6);		// add vertex 6 at time 6

	vertices.at(0) = 7;
	tc.add_insertion(vertices, 7);		// add vertex 7 at time 7

	vertices.at(0) = 8;
	tc.add_insertion(vertices, 8);		// add vertex 8 at time 8

	vertices.at(0) = 9;
	tc.add_insertion(vertices, 9);		// add vertex 9 at time 9

	tc.add_contraction(5, 0, 10);		// add contraction of 5 and 0 at time 10, 5 is not valid from now on.

	vertices.at(0) = 10;
	tc.add_insertion(vertices, 11);		// add vertex 10 at time 11

	vertices.at(0) = 0;
	vertices.push_back(7);
	tc.add_insertion(vertices, 12);		// add edge (0,7) at time 12

	vertices.at(0) = 0;
	vertices.at(1) = 8;
	tc.add_insertion(vertices, 13);		// add edge (0,8) at time 13

	vertices.at(0) = 1;
	vertices.at(1) = 3;
	tc.add_insertion(vertices, 14);		// add edge (1,3) at time 14

	tc.add_contraction(6, 1, 15);		// add contraction of 6 and 1 at time 15, 6 is not valid from now on.

	vertices.pop_back();
	vertices.at(0) = 11;
	tc.add_insertion(vertices, 16);		// add vertex 11 at time 16

	vertices.at(0) = 1;
	vertices.push_back(7);
	tc.add_insertion(vertices, 17);		// add edge (1,7) at time 17

	tc.add_contraction(8, 1, 18);		// add contraction of 8 and 1 at time 18, 8 is not valid from now on.

	// ...

	outputFile->close();
}

void write_tc_results_to_file_as_vertices(Hash_complex *complex, std::vector<Hash_complex::vertex> &simplex, double timestamp){
	std::vector<Hash_complex::vertex>::size_type size = simplex.size();
	*outputFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (size - 1) << " ";
	for (std::vector<Hash_complex::vertex>::size_type i = 0; i < size; i++){
		*outputFile << simplex.at(i) << " ";
	}
	*outputFile << timestamp << std::endl;
}

void write_tc_results_to_file_as_faces(Hash_complex *complex, std::vector<Hash_complex::vertex> &simplex, double timestamp){
	std::vector<Hash_complex::vertex>::size_type size = simplex.size();
	*outputFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << (size - 1) << " ";
	if (size > 1){
		std::vector<Hash_complex::index> boundary;
		complex->get_boundary(simplex, &boundary);
		for (std::vector<Hash_complex::index>::size_type i = 0; i < size; i++){
			*outputFile << boundary.at(i) << " ";
		}
	}
	*outputFile << timestamp << std::endl;
}

