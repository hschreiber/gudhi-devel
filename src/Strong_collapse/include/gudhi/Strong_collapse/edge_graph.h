/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STRONG_COLLAPSE_EDGE_GRAPH_H_
#define STRONG_COLLAPSE_EDGE_GRAPH_H_

#include <algorithm>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include <gudhi/graph_simplicial_complex.h>

namespace Gudhi {

namespace strong_collapse {

template<typename SimplicialComplexForProximityGraph>
class Edge_graph {
public:
	using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;
	using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
	using Filtered_edge = std::tuple<Filtration_value, Vertex_handle, Vertex_handle>;
	using Filtered_edge_set = std::vector<Filtered_edge>;

	Edge_graph(Proximity_graph<SimplicialComplexForProximityGraph>& graph);
	Edge_graph(const Edge_graph& graph);

	std::size_t size() const;
	Filtration_value get_filtration_min() const;
	Filtration_value get_filtration_max() const;
	Filtered_edge_set sub_filter_edges_by_filtration(Filtration_value new_threshold);
	Filtered_edge_set sub_filter_edges_by_index(std::size_t new_index);
	Filtration_value get_filtration_at(std::size_t index) const;

private:
	Filtered_edge_set edges_;
};

template<typename SimplicialComplexForProximityGraph>
inline Edge_graph<SimplicialComplexForProximityGraph>::Edge_graph(Proximity_graph<SimplicialComplexForProximityGraph> &graph)
{
	using edge_iterator = typename boost::graph_traits<Proximity_graph<SimplicialComplexForProximityGraph> >::edge_iterator;

	std::pair<edge_iterator, edge_iterator> itPair = boost::edges(graph);

	for (edge_iterator edge_it = itPair.first; edge_it != itPair.second; ++edge_it) {
		edges_.push_back(Filtered_edge(
							 boost::get(edge_filtration_t(), graph, *edge_it),
							 boost::source(*edge_it, graph),
							 boost::target(*edge_it, graph)
							 ));
	}

	std::sort(edges_.begin(),
			  edges_.end(),
			  [](const Filtered_edge& e1, const Filtered_edge& e2){
				  return std::get<0>(e1) < std::get<0>(e2);
			  }
	);
}

template<typename SimplicialComplexForProximityGraph>
inline Edge_graph<SimplicialComplexForProximityGraph>::Edge_graph(const Edge_graph &graph)
	: edges_(graph.edges_)
{}

template<typename SimplicialComplexForProximityGraph>
inline std::size_t Edge_graph<SimplicialComplexForProximityGraph>::size() const
{
	return edges_.size();
}

template<typename SimplicialComplexForProximityGraph>
inline typename Edge_graph<SimplicialComplexForProximityGraph>::Filtration_value
Edge_graph<SimplicialComplexForProximityGraph>::get_filtration_min() const
{
	return std::get<0>(edges_.front());
}

template<typename SimplicialComplexForProximityGraph>
inline typename Edge_graph<SimplicialComplexForProximityGraph>::Filtration_value
Edge_graph<SimplicialComplexForProximityGraph>::get_filtration_max() const
{
	return std::get<0>(edges_.back());
}

template<typename SimplicialComplexForProximityGraph>
inline typename Edge_graph<SimplicialComplexForProximityGraph>::Filtered_edge_set
Edge_graph<SimplicialComplexForProximityGraph>::sub_filter_edges_by_filtration(Filtration_value new_threshold)
{
	if (new_threshold > get_filtration_max()) return edges_;

	Filtered_edge_set res;

	auto it = edges_.begin();
	while (std::get<0>(*it) < new_threshold) {
		res.push_back(*it);
		++it;
	}

#ifdef DEBUG_TRACES
	if (!res.empty()) {
		Filtered_edge back = res.back();
		std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration threshold = "
				<< std::get<0>(back) << " - size = " << res.size() << " - back = "
				<< std::get<0>(back) << " - u = " << std::get<1>(back)
				<< " - v = " << std::get<2>(back) << std::endl;
	} else {
		std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration is empty " << std::endl;
	}
#endif  // DEBUG_TRACES

	return res;
}

template<typename SimplicialComplexForProximityGraph>
inline typename Edge_graph<SimplicialComplexForProximityGraph>::Filtered_edge_set
Edge_graph<SimplicialComplexForProximityGraph>::sub_filter_edges_by_index(std::size_t new_index)
{
#ifdef DEBUG_TRACES
	if (new_index < edges_.size()) {
		Filtered_edge_set output(edges_.begin(), edges_.begin() + new_index + 1);
		Filtered_edge back(output.back());
		std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration threshold = "
				<< std::get<0>(back) << " - size = " << output.size() << " - back = "
				<< std::get<0>(back) << " - u = " << std::get<1>(back)
				<< " - v = " << std::get<2>(back) << std::endl;
	} else {
		std::cout << "Filtered_edges_vector::sub_filter_edges_by_index is empty " << std::endl;
	}
#endif  // DEBUG_TRACES

	if (new_index >= edges_.size()) return edges_;

	return Filtered_edge_set(edges_.begin(), edges_.begin() + new_index + 1);
}

template<typename SimplicialComplexForProximityGraph>
inline typename Edge_graph<SimplicialComplexForProximityGraph>::Filtration_value
Edge_graph<SimplicialComplexForProximityGraph>::get_filtration_at(std::size_t index) const
{
	return std::get<0>(edges_[index]);
}

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // STRONG_COLLAPSE_EDGE_GRAPH_H_
