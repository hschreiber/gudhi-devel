/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "tower_assembler"
#include <boost/test/unit_test.hpp>

#include <boost/test/test_tools.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <tuple>
#include <cstdio>
#include <limits>
#include <cmath>

#include <gudhi/Strong_collapse/Flag_complex_tower_assembler.h>
#include <gudhi/Strong_collapse/Flag_complex_sparse_matrix.h>
#include <gudhi/Unitary_tests_utils.h>

// Type definitions
struct Simplicial_complex {
	using Vertex_handle = int;
	using Filtration_value = double;
};

using Flag_complex_sparse_matrix = Gudhi::strong_collapse::Flag_complex_sparse_matrix<Simplicial_complex>;
using Filtered_sorted_edge_list = Flag_complex_sparse_matrix::Filtered_sorted_edge_list;
using Edge_list = Flag_complex_sparse_matrix::Edge_list;
using Reduction_map = Flag_complex_sparse_matrix::Reduction_map;
using Flag_complex_tower_assembler = Gudhi::strong_collapse::Flag_complex_tower_assembler<Simplicial_complex>;

#define GUDHI_TEST_DISTANCE_MATRIX_CLOSE(first_collection, second_collection) { \
	BOOST_CHECK(first_collection.size() == second_collection.size()); \
	auto second_iter = std::begin(second_collection); \
	int i=0, j=0; \
	for (auto &line : first_collection) { \
	auto second_iter_iter = std::begin(*second_iter); \
	for (auto value : line) { \
	BOOST_CHECK_MESSAGE(std::fabs(value - *second_iter_iter) <= std::numeric_limits<double>::epsilon(), \
	"Failed for value [" << i << "][" << j << "] - " << value << " versus " << *second_iter_iter); \
	second_iter_iter++; \
	j++; \
	} \
	second_iter++; \
	i++; j=0; \
	} \
	}

BOOST_AUTO_TEST_CASE(tower_assembler_strong_collapse) {
	// Clean produced files
	remove("./tower_assembler_1.txt");
	remove("./tower_assembler_2.txt");
	remove("./tower_assembler_3.txt");
	remove("./tower_assembler_4.txt");
	remove("./tower_assembler_5.txt");

	Filtered_sorted_edge_list input_edges;
	input_edges.push_back(std::make_tuple(1., 1, 2));
	input_edges.push_back(std::make_tuple(1., 2, 3));
	input_edges.push_back(std::make_tuple(1., 3, 4));
	input_edges.push_back(std::make_tuple(1., 4, 1));

	// Edges:
	//  1 o---o 2
	//    |   |
	//    |   |
	//    |   |
	//  4 o---o 3
	// No collapse is performed
	Flag_complex_sparse_matrix mat_coll_1(4, input_edges);

	mat_coll_1.strong_collapse();

	input_edges.push_back(std::make_tuple(1.5, 1, 3));
	input_edges.push_back(std::make_tuple(1.5, 2, 4));

	// Edges:
	//  1 o---o 2
	//    |\ /|
	//    | X |
	//    |/ \|
	//  4 o---o 3
	// All is collapsed on 1
	Flag_complex_sparse_matrix mat_coll_2(4, input_edges);

	mat_coll_2.strong_collapse();

	Flag_complex_tower_assembler tower_assembler_1(4);

	tower_assembler_1.build_tower_for_two_complexes(mat_coll_1,
													mat_coll_2,
													10.,
													"./tower_assembler_1.txt");
	Flag_complex_sparse_matrix::Distance_matrix sparse_distances_1 = tower_assembler_1.distance_matrix();

	std::cout << "#1. sparse_distances size =" << sparse_distances_1.size() << std::endl;
	for (auto &line : sparse_distances_1) {
		for (auto value : line) {
			std::cout << value << ", ";
		}
		std::cout << std::endl;
	}

	input_edges.push_back(std::make_tuple(2., 2, 5));
	input_edges.push_back(std::make_tuple(2., 5, 6));
	input_edges.push_back(std::make_tuple(2., 3, 6));

	// Edges:
	//        2
	//  1 o---o---o 5
	//    |\ /|   |
	//    | X |   |
	//    |/ \|   |
	//  4 o---o---o 6
	//        3
	// Check 4 and 1 are collapsed on 2
	Flag_complex_sparse_matrix mat_coll_3(6, input_edges);

	mat_coll_3.strong_collapse();

	Flag_complex_tower_assembler tower_assembler_2(6);

	tower_assembler_2.build_tower_for_two_complexes(mat_coll_2,
													mat_coll_3,
													10.,
													"./tower_assembler_2.txt");
	Flag_complex_sparse_matrix::Distance_matrix sparse_distances_2 = tower_assembler_2.distance_matrix();

	std::cout << "#2. sparse_distances size =" << sparse_distances_2.size() << std::endl;
	for (auto &line : sparse_distances_2) {
		for (auto value : line) {
			std::cout << value << ", ";
		}
		std::cout << std::endl;
	}

	input_edges.push_back(std::make_tuple(2.5, 2, 6));
	input_edges.push_back(std::make_tuple(2.5, 5, 3));

	// Edges:
	//        2
	//  1 o---o---o 5
	//    |\ /|\ /|
	//    | X | X |
	//    |/ \|/ \|
	//  4 o---o---o 6
	//        3
	// all is collapsed on 2
	Flag_complex_sparse_matrix mat_coll_4(6, input_edges);

	mat_coll_4.strong_collapse();

	Flag_complex_tower_assembler tower_assembler_3(6);

	tower_assembler_3.build_tower_for_two_complexes(mat_coll_3,
													mat_coll_4,
													10.,
													"./tower_assembler_3.txt");
	Flag_complex_sparse_matrix::Distance_matrix sparse_distances_3 = tower_assembler_3.distance_matrix();

	std::cout << "#3. sparse_distances size =" << sparse_distances_3.size() << std::endl;
	for (auto &line : sparse_distances_3) {
		for (auto value : line) {
			std::cout << value << ", ";
		}
		std::cout << std::endl;
	}
	GUDHI_TEST_DISTANCE_MATRIX_CLOSE(sparse_distances_1, sparse_distances_3);

	input_edges.push_back(std::make_tuple(3., 3, 7));
	input_edges.push_back(std::make_tuple(3., 7, 8));
	input_edges.push_back(std::make_tuple(3., 4, 8));

	// Edges:
	//        2
	//  1 o---o---o 5
	//    |\ /|\ /|
	//    | X | X |
	//    |/ \|/ \|
	//  4 o---o---o 6
	//    |   | 3
	//    |   |
	//    |   |
	//  8 o---o 7
	// 1, 2, 5 and 6 are collapsed on 3
	Flag_complex_sparse_matrix mat_coll_5(8, input_edges);

	mat_coll_5.strong_collapse();

	Flag_complex_tower_assembler tower_assembler_4(8);

	tower_assembler_4.build_tower_for_two_complexes(mat_coll_4,
													mat_coll_5,
													10.,
													"./tower_assembler_4.txt");
	Flag_complex_sparse_matrix::Distance_matrix sparse_distances_4 = tower_assembler_4.distance_matrix();

	std::cout << "#4. sparse_distances size =" << sparse_distances_4.size() << std::endl;
	for (auto &line : sparse_distances_4) {
		for (auto value : line) {
			std::cout << value << ", ";
		}
		std::cout << std::endl;
	}
	GUDHI_TEST_DISTANCE_MATRIX_CLOSE(sparse_distances_2, sparse_distances_4);

	input_edges.push_back(std::make_tuple(3., 3, 8));
	input_edges.push_back(std::make_tuple(3., 4, 7));

	// Edges:
	//        2
	//  1 o---o---o 5
	//    |\ /|\ /|
	//    | X | X |
	//    |/ \|/ \|
	//  4 o---o---o 6
	//    |\ /| 3
	//    | X |
	//    |/ \|
	//  8 o---o 7
	// all is collapsed on 3
	Flag_complex_sparse_matrix mat_coll_6(8, input_edges);

	mat_coll_6.strong_collapse();

	Flag_complex_tower_assembler tower_assembler_5(8);

	tower_assembler_5.build_tower_for_two_complexes(mat_coll_5,
													mat_coll_6,
													10.,
													"./tower_assembler_5.txt");
	Flag_complex_sparse_matrix::Distance_matrix sparse_distances_5 = tower_assembler_5.distance_matrix();

	std::cout << "#5. sparse_distances size =" << sparse_distances_5.size() << std::endl;
	for (auto &line : sparse_distances_5) {
		for (auto value : line) {
			std::cout << value << ", ";
		}
		std::cout << std::endl;
	}
	GUDHI_TEST_DISTANCE_MATRIX_CLOSE(sparse_distances_3, sparse_distances_5);
}

