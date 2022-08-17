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
	std::cout << "  ./example_tower_from_file_write_filtration_into_file input_file_name output_file_name" << std::endl;
}

std::ofstream *outputFile;		/**< Pointer to output file. */

void write_tc_results_to_file_as_vertices(Hash_complex *complex, std::vector<Hash_complex::vertex> &simplex, double timestamp);	/**< Example 1 of callback function for Tower_converter. */
void write_tc_results_to_file_as_faces(Hash_complex *complex, std::vector<Hash_complex::vertex> &simplex, double timestamp);	/**< Example 2 of callback function for Tower_converter. */

int main(int argc, char *argv[])
{
    if (argc != 3){
		print_usage();
		return 0;
    }

    std::ifstream file(argv[1]);
	outputFile = new std::ofstream(argv[2]);

	if (!outputFile->is_open()){
		std::cout << "Unable to open output file." << std::endl;
		if (!file.is_open()){
			std::cout << "Unable to open input file." << std::endl;
			file.setstate(std::ios::failbit);
		}
		return 0;
	}
	if (!file.is_open()){
		std::cout << "Unable to open input file." << std::endl;
		file.setstate(std::ios::failbit);
		outputFile->close();
		return 0;
	}

	Tower_converter<Hash_complex> tc(write_tc_results_to_file_as_vertices);	// output callback function ; another example: write_tc_results_to_file_as_faces.
	file >> tc;		// >> function in gudhi/tc_reading_utilities.h, see documentation for input file format.

	file.close();
	outputFile->close();
	return 0;
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



