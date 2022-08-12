#include <gudhi/Flag_complex_strong_collapse.h>  // for Flag_complex_strong_collapse
#include <gudhi/graph_simplicial_complex.h>  // for Filtered_edges_vector
#include <gudhi/distance_functions.h>  // for Euclidean_distance
#include <gudhi/Strong_collapse/edge_graph.h>

#include <iostream>
#include <vector>
#include <iomanip>  // for std::setw, std::setprecision

struct Simplicial_complex {
	using Vertex_handle = short;
	using Filtration_value = double;
};

using Point = std::vector<Simplicial_complex::Filtration_value>;
using Point_cloud = std::vector<Point>;
using Proximity_graph = Gudhi::Proximity_graph<Simplicial_complex>;
using Edge_graph = Gudhi::strong_collapse::Edge_graph<Simplicial_complex>;
using Flag_complex_strong_collapse = Gudhi::strong_collapse::Flag_complex_strong_collapse<Simplicial_complex>;

int main(int argc, char **argv) {
	// ----------------------------------------------------------------------------
	// Init of a list of points from a small molecule
	// ----------------------------------------------------------------------------
	Point_cloud points;
	points.push_back({1.0, 1.0});
	points.push_back({7.0, 0.0});
	points.push_back({4.0, 6.0});
	points.push_back({9.0, 6.0});
	points.push_back({0.0, 14.0});
	points.push_back({2.0, 19.0});
	points.push_back({9.0, 17.0});

	// ----------------------------------------------------------------------------
	// Edge graph construction from the list of points
	// ----------------------------------------------------------------------------
	Proximity_graph prox_graph =
			Gudhi::compute_proximity_graph<Simplicial_complex>(
				points,
				12.,  // means all edges
				Gudhi::Euclidean_distance());

	Edge_graph edge_graph(prox_graph);

	std::cout << "Edge graph contains " << edge_graph.size() << " edges. ";
	std::cout << "Minimal edge length is '" << edge_graph.get_filtration_min() << "'. ";
	std::cout << "Maximal edge length is '" << edge_graph.get_filtration_max() << "'." << std::endl;

	Flag_complex_strong_collapse flag_complex_exact_version(points.size(), prox_graph);
	std::cout << "Distance matrix after strong collapse computation." << std::endl;
	auto distance_matrix = flag_complex_exact_version.get_distance_matrix();

	for (auto &line : distance_matrix) {
		for (auto value : line) {
			std::cout << std::setw(12) << std::setprecision(4) << value;
		}
		std::cout << std::endl;
	}
	return 0;
}
